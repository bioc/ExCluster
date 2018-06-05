ExCluster <- function(exonCounts=NULL, cond.Nums=NULL, annot.GFF = NULL, GFF.File=NULL, out.Dir=NULL, result.Filename=NULL,CombineExons=FALSE,plot.Results=FALSE,FDR.cutoff=0.05){

    ### ExCluster progression message
    cat("","Initializing ...",sep="\n")

    ### check to make sure a normalized exon count matrix was provided as input
    if (is.null(exonCounts) == TRUE){
        stop(call = "You did provide a normalized exon count matrix to the exonCounts argument.Please obtain normalized exon counts from processCounts() in the form of a global variable, and assign said variable to exonCounts.
             For example, if you followed the vignette, you should have a normCounts variable. Please specify the exonCounts=normCounts argument when running this function.")
    }

    ### ensure the number of columns in the exonCounts parameter equal the number of cond.Nums
    if (ncol(exonCounts) != length(cond.Nums)){
        stop(call = "Your number of exonCounts columns did not equal the length of your number of condition identifiers. Please ensure that you entered correct condition identifiers to cond.Nums, such as cond.Nums=c(1,1,1,2,2,2) for 3 repliates per 2 conditions.
    Additionally, please also ensure your normalized counts provided to exonCounts have only count data columns, and that Gene:ExonBin identifiers are provided as rownames.
    If you are unsure, please refer to the vignette, or re-run processCounts and assign its output to a variable, and then assign that variable to the exonCounts argument.")
    }

    ### check to ensure that there are no duplicated data columns in exonCounts
    CheckDups <- duplicated(t(exonCounts))
    if (any(CheckDups) == TRUE){
        DupCols <- which(CheckDups == TRUE)
        print(paste("Columns",DupCols,"were duplicated in your normalized exon count data. Please re-check your BAM file names and re-run processCounts."))
        stop(call = "One of your normalized exon count columns was duplicated! If this is a mistake, please ensure you specified the correct BAM files to processCounts, and re-run that function.
    However, if you are attempting to duplicate count columns to bypass ExClusters requirement for biological replicates, please be aware that this will not work well. ExCluster computes variances per condition, and one of your conditions will have zero variance. This will results in high false positive rates. ")
    }

    ### check to make sure either a GFF data structure or GFF file path were specified
    if (is.null(annot.GFF) == TRUE){
        if (is.null(GFF.File) == TRUE){
            stop(call = "You must provide either GFF annotation data or a GFF file path. Please assign the annot.GFF argument a variable name containing the GFF data, such as annot.GFF=GFF. Alternatively, please specify a full GFF file path, such as: GFF.File=/Users/username/path/to/file.gff")
        }else{
            ### Test if GFF.File works
            if(file.exists(GFF.File) != TRUE){
                stop(call = "The GFF file path you entered did not exist. Please verify that you have assigned the exactly correct path to the GFF.File argument.
    For example, your argument when running ExCluster should look something like this: GFF.File=/Users/username/path/to/file.gff")
            }
        }
    }

    #####################################################################################################
    #######################################  INPUT PARAMETERS  ##########################################
    #####################################################################################################

    ## clustering method = euclidean
    DISTANCE <- "euclidean"

    ## linkage == complete
    LINKAGE <- "complete"


    #####################################################################################################
    ############################################  FUNCTIONS  ############################################
    #####################################################################################################


    #### Geometric mean function
    gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }

    ### function to determine number of clusters
    Determine_HC_Distance <- function(x){
        which(HC_cutree%in%x)
    }

    ### function to compute the distances between all  clusters (apply max to its output)
    Compare_HC_Distances <- function(z){
        ## convert ES_mat matrix into symmetric dist matrix
        ES_dist_mat <- as.matrix(dist(ES_mat,method=DISTANCE))
        lapply(1:NumClusters,FUN=function(y){
            TempMat <- ES_dist_mat[unlist(HC_indices[[z]]),unlist(HC_indices[[y]])]
            TempMat <- unlist(lapply(TempMat, function(w) ifelse(w==0, 0.0000001, w)))
            gm_mean(TempMat)})
    }

    #Compute effect size
    Calc.effect.size <- function(Mean1,Mean2,Var1,Var2){
        Effect_size <- (Mean1-Mean2)/(sqrt((Var1+Var2)/2))
    }

    #function to generate all pairwise effect size distances
    Generate_ES_dists <- function(x){
        Calc.effect.size(Mean1[x],Mean2[x],Var1[x],Var2[x])
    }

    Permuted_ES_nullhypo <- function(r, Sim_log2FC=NULL, Sim_log2var=NULL, NumRows=NULL){
        Mean1 <- rep(Sim_log2FC[r,], NumRows)
        Mean2 <- rep(Sim_log2FC[r,], each=NumRows)
        Var1 <- rep(Sim_log2var[r,], NumRows)
        Var2 <- rep(Sim_log2var[r,], each=NumRows)
        ES_permuted <- matrix(unlist(lapply(1:(NumRows^2),
            FUN=function(x){(Mean1[x]-Mean2[x])/(sqrt((Var1[x]+Var2[x])/2))})),
            nrow=NumRows,byrow=TRUE)
        #ES_permuted <- matrix(Generate_ES_dists(1:(NumRows^2)),nrow=NumRows,byrow=TRUE)
        if (NumRows > 3){
            ##Generate hclust object
            HC <- hclust(dist(ES_permuted,DISTANCE),method=LINKAGE)
            ##cut tree
            HC_cutree <- cutree(HC,k=NumClusters)
            ## find row indices per cluster
            HC_indices <- lapply(1:NumClusters,FUN=function(y){which(HC_cutree%in%y)})
        }else{
            HC_indices <- lapply(1:NumRows,FUN=function(y){y})
        }
        ### generate distance matrix from distance object
        ES_SimMat <- as.matrix(dist(ES_permuted,method=DISTANCE))
        ES_Clust_Dists <- matrix(unlist(lapply(1:NumClusters,FUN=function(z){
            lapply(1:NumClusters,FUN=function(y){
                TempMat <- ES_SimMat[unlist(HC_indices[[z]]),unlist(HC_indices[[y]])]
                TempMat <- unlist(lapply(TempMat, function(w) ifelse(w==0, 0.0000001, w)))
                gm_mean(TempMat)
            })
        })),nrow=NumClusters,byrow=TRUE)
        return(max(apply(ES_Clust_Dists,1,sum)/(NumClusters-1)))
    }

    ### Function for running null hypothesis non-parametric iterations
    generate_Nullhypo_Dists <- function(NumIterations=NULL,gene_FC=NULL,Gene_R1=NULL,Gene_R2=NULL, NumRows=NULL){

        ## rows = genes, columns = exons
        Sim_log2FC <- matrix(NA,nrow=NumIterations,ncol=NumRows)
        Sim_log2var <- matrix(NA,nrow=NumIterations,ncol=NumRows)

        for (l in 1:NumRows){
            SimulatedMatrix1 <- matrix(abs(SimulatedMatrix1 <- rnorm(n = NumIterations*length(log2Indices1), mean=log2((apply(log2Clusters[rows[l],c(log2Indices1,log2Indices2)],1,mean)*2*Gene_R1)+1), sd=sqrt(log2Clusters$var1[rows[l]]))),nrow=NumIterations)
            SimulatedMatrix2 <- matrix(abs(SimulatedMatrix2 <- rnorm(n = NumIterations*length(log2Indices2), mean=log2((apply(log2Clusters[rows[l],c(log2Indices1,log2Indices2)],1,mean)*2*Gene_R2)+1), sd=sqrt(log2Clusters$var2[rows[l]]))),nrow=NumIterations)

            Sim_log2FC[,l] <- apply(SimulatedMatrix2,1,mean) - apply(SimulatedMatrix1,1,mean)
            Sim_log2var[,l] <- apply(SimulatedMatrix2,1,var) + apply(SimulatedMatrix1,1,var)
        }

        Nonparam_Clust_Dists <- unlist(lapply(1:NumIterations,Permuted_ES_nullhypo,Sim_log2FC=Sim_log2FC, Sim_log2var=Sim_log2var, NumRows=NumRows))
        return(Nonparam_Clust_Dists)
    }

    ### Note this function was adapted in full or in part from the NbClust R packaage
    # NbClust offers no warranty on this function as per their GPL-2 license
    centers <- function(cl, x) {
        x <- as.matrix(x)
        n <- length(cl)
        k <- max(cl)
        centers <- matrix(nrow = k, ncol = ncol(x))
        {
            for (i in 1:k) {
                for (j in 1:ncol(x)) {
                    centers[i, j] <- mean(x[cl == i, j])
                }
            }
        }
        return(centers)
    }

    ### Note this function was adapted in full or in part from the NbClust R packaage
    # NbClust offers no warranty on this function as per their GPL-2 license
    ### SD index scattering computation
    SD.Average.scattering <- function(cl, x) {
        x <- as.matrix(x)
        n <- length(cl)
        k <- max(cl)
        centers.matrix <- centers(cl, x)
        cluster.size <- numeric(0)
        variance.clusters <- matrix(0, ncol = ncol(x), nrow = k)
        var <- matrix(0, ncol = ncol(x), nrow = k)
        for (u in 1:k) cluster.size[u] <- sum(cl == u)
        for (u in 1:k) {
            for (j in 1:ncol(x)) {
                for (i in 1:n) {
                    if (cl[i] == u)
                        variance.clusters[u, j] <- variance.clusters[u, j] + (x[i, j] - centers.matrix[u, j])^2
                }
            }
        }
        for (u in 1:k) {
            for (j in 1:ncol(x)) variance.clusters[u, j] = variance.clusters[u, j]/cluster.size[u]
        }
        variance.matrix <- numeric(0)
        for (j in 1:ncol(x)) variance.matrix[j] = var(x[, j]) * (n - 1)/n
        Somme.variance.clusters <- 0
        for (u in 1:k) Somme.variance.clusters <- Somme.variance.clusters + sqrt((variance.clusters[u, ] %*% (variance.clusters[u, ])))
        stdev <- (1/k) * sqrt(Somme.variance.clusters)
        scat <- (1/k) * (Somme.variance.clusters/sqrt(variance.matrix %*% variance.matrix))
        scat <- list(stdev = stdev, centers = centers.matrix, variance.intraclusters = variance.clusters, scatt = scat)
        return(scat)
    }

    ### Note this function was adapted in full or in part from the NbClust R packaage
    # NbClust offers no warranty on this function as per their GPL-2 license
    #### SD Index distance computation
    SD.Distance <- function(cl,x){
        Clust_Dists <- array()
        Dist_Counter <- 1
        MAX_cl <- max(cl)
        for (m in 1:MAX_cl){
            N <- 1:MAX_cl
            N <- N[-c(m)]
            for (n in N){
                Clust_Dists[Dist_Counter] <- mean(unlist(x[which(cl%in%m),which(cl%in%n)]))
                Dist_Counter <- Dist_Counter + 1
            }
        }
        Mean_Dist <- as.numeric(mean(Clust_Dists))
        return(Mean_Dist)
    }

    ### Note this function was adapted in full or in part from the NbClust R packaage
    # NbClust offers no warranty on this function as per their GPL-2 license
    #### SD index final computation, combining scattering + distance results
    SD.Index2 <- function(x, clmin, cl) {
        x <- as.matrix(x)
        Scatt <- SD.Average.scattering(cl, x)$scatt
        Dis0 <- SD.Distance(cl, x)
        Alpha <- SD.Distance(clmin, x)
        SD.indices2 <- Dis0/Alpha + (1-Scatt)
        return(SD.indices2)
    }

    ### Note this function was adapted in full or in part from the NbClust R packaage
    # NbClust offers no warranty on this function as per their GPL-2 license
    Cutree_SD.Index <- function(){
        ## hard coded at minimum of 2 clusters
        SD_values <- array()
        for(C in 2:NRows){
            ## cut tree based on observed number of clusters C
            cl1 <- cutree(HC,k=C)
            ## cut tree into the minimum number of clusters (k=2)
            clmin <- cutree(HC, k=2)
            x <- as.matrix(dist(ES_mat,method = DISTANCE))
            SD_values[(C-1)] <- as.numeric(SD.Index2(x,clmin,cl1))
            if (SD_values[(C-1)] < max(SD_values)){
                break
            }
        }
        ### estimated number clusters is the maximized SD_value, + 1 (because we start at 2 clusters)
        Estimated_Clustnum <- (C-1)
        res <- cutree(HC,k=Estimated_Clustnum)
        return(res)
    }

    ### function for estimating the lowest null hypothesis p-value (very rough estimate)
    # x = p-value vector
    EstNullHypo <- function(x){
        # number of p-values in every interval of 0.02
        NumEstNull <- array()
        for (i in 6:56){
            Lower <- (i)/400
            Upper <- (i+16)/400
            NumEstNull[i+1] <- length(x[which(x >= Lower & x < Upper)])
        }
        # check to see if p-values have hit a low plateau (null hypothesis)
        for (i in 6:40){
            # Check maximum number of pvals in the 4 below the current i
            NumBelow <- NumEstNull[(i+1):(i+5)]
            # Convert these to TRUE/FALSE
            NumBelowBoolean <- NumEstNull[i] > NumBelow
            # Now see how many of these bins the current NumEstNull[i] is less than
            LessThan <- length(which(NumBelowBoolean == FALSE))
            if (LessThan >= 2) break
        }
        # now your i value divided by 400 is your pvalue cutoff
        PValCutoff <-i/400
        return(PValCutoff)
    }

    ## function for correcting FDR -- to be run with lapply()
    # x = the p-value vector
    # y = the number of estimated True p-values
    # z = the number of estimated Null Hypothesis pvalues
    FDRcalc <- function(x,y,z){
        # set p FDR values array to be filled
        FDRvalues <- NULL
        # compute the prior odds ratio (NumTrue/NumNull)
        PriorOdds <- y/z

        # now loop through pvalues
        for (n in 1:length(x)){
            # check to make sure that we aren't on the first value, and the previous value isn't FDR == 1
            if (n > 1 && FDRvalues[(n-1)] >= 0.99) {
                    FDRvalues[n] <- 1
            }else{
                ### now generate a null hypothesis 2 times to estimate the NumExp and Num Obs
                NumExp <- NULL
                NumObs <- NULL
                for (h in 1:3){
                    # generate pure null hypothesis distribution
                    NullHypo <- sort(runif(z,0,1))
                    # number of expected p-values based on the null hypothesis distribution
                    NumExp[h] <- floor(length(which(NullHypo <= x[n])))
                    # number of observed TRUE p-values less than or equal to x[n]
                    NumObs[h] <- length(which(x <= x[n])) - NumExp[h]
                    # make sure NumObs is not negative
                    NumObs[h] <- max(NumObs[h],0)
                }
                # now take the average of these estimates
                NumExp <- NumExp[1]
                NumObs <- NumObs[1]

                # compute the likelihood ratio
                L_Ratio <- NumObs/NumExp
                # compute posterior odds ratio, which is the PriorOdds * the likelihood ratio
                PostOdds <- PriorOdds * L_Ratio
                # compute posterior probability
                PostProb <- 1/(1+(1/PostOdds))
                # now compute FDR value, which is 1 minus the posterior probability
                FDRvalues[n] <- 1 - PostProb
            }
        }

        ### now loop back through FDRs and make sure that no FDR is lower than higher ranked values, stop when FDRs==1
        # find the limit at which FDR==1
        FDRlim <- min(which(FDRvalues == 1),(length(x)-2))
        # loop through each FDR value
        for (n in 2:FDRlim[1]){
            # make sure equal p-values have the same FDR for pvals below
            if (x[n] == x[(n-1)]){
                if (FDRvalues[(n-1)] > FDRvalues[n]){
                    FDRvalues[(n-1)] <- FDRvalues[n]
                }
            }
            # make sure equal p-values have the same FDR for pvals above
            if (x[n] == x[(n+1)]){
                if (FDRvalues[n] > FDRvalues[(n+1)]){
                    FDRvalues[n] <- FDRvalues[(n+1)]
                }
            }

            # compute the minimum FDR above FDRvalues[n]
            MaxN <- min((n+500),FDRlim[1])
            MinFDR <- min(FDRvalues[(n-1):MaxN])
            # now replace the FDRvalues[n] if FDRvalues[n] > Min FDR
            if (FDRvalues[(n-1)] > MinFDR){
                FDRvalues[(n-1)] <- MinFDR
            }
        }
        return(FDRvalues)
    }


    #####################################################################################################
    #################################  READ COUNT & DATA PREPARATION  ###################################
    #####################################################################################################

    ### ExCluster progression message
    cat('','Reading in data ...','',sep="\n")

    log2FC <- exonCounts
    rm(exonCounts)

    ## grabbing unique condition numbers (should only ever be two conditions)
    uCondNums <- unique(cond.Nums)

    ## ensure exactly 2 conditions
    if (length(uCondNums) != 2){
        stop(call = "You did not specify exactly 2 condition numbers to the cond.Nums argument. For example, if your data contains two conditions with three replicates each, in order, please use: cond.Nums=c(1,1,1,2,2,2)
             As another example, if your count data alternates columns with cond1, cond2, cond1, cond2, etc., please use: cond.Nums=c(1,2,1,2,1,2)")
    }

    ## condition indices
    Indices1 <- grep(uCondNums[1],cond.Nums)
    Indices2 <- grep(uCondNums[2],cond.Nums)

    ## Ensure that both conditions have at least 2 replicates provided, or ExCluster will not function.
    if (length(Indices1) == 1){
        stop(call = "You only provided biological replicate for your first condition.
    ExCluster, unfortunately, requires at least 2 biological replicates per condition to function.
    If this was a mistake, please re-enter your cond.Nums argument with at least 2 replicates per condition number.")
    }
    if (length(Indices2) == 1){
        stop(call = "You only provided biological replicate for your second condition.
    ExCluster, unfortunately, requires at least 2 biological replicates per condition to function.
    If this was a mistake, please re-enter your cond.Nums argument with at least 2 replicates per condition number.")
    }

    ## Grab number of replicates from Cond1 and Cond2
    NumReps1 <- length(Indices1)
    NumReps2 <- length(Indices2)

    ######################### READ IN GFF FILE IF NECESSARY #########################

    ### Read in the GFF file if the annot.GFF argument was not specified
    if (is.null(annot.GFF) == TRUE){
        if (is.null(GFF.File) == FALSE){
            ### Read in GFF file and grab exon bins (THIS MUST BE THE EXACT GFF FILE USED TO GENERATE READ COUNTS!!!)
            annot.GFF <- read.table(file=GFF.File,sep="\t",stringsAsFactors=FALSE)
        }
    }else{
        if (ncol(annot.GFF) != 10){
            ### inform the user that the GFF file was not properly formatted
            stop(call="You did not enter a valid GFF annotation data frame or a proper path to your GFF file.
    It is possible that your GFF file is not properly formatted with the expected 10 columns.
    Please re-run GFF_convert and then re-assign either annot.GFF or GFF.File when calling ExCluster. See the ExCluster manual and vignette for more information.")
        }
    }

    ######################### NOW BEAT DATA INTO SHAPE FOR log2FC and annot.GFF #########################

    # first sort log2FC on EnsG:exon_bin
    log2FC <- log2FC[order(rownames(log2FC)),]
    # do the same for annot.GFF
    annot.GFF <- annot.GFF[order(annot.GFF$V2),]

    # now copy log2FC to newlog2FC, so we can work with newlog2FC and preserve log2FC
    newlog2FC <- log2FC

    ## now assign Gene names and exon numbers for later joining to log2FC -- used for plotting exons
    log2Bins <- gsub(pattern = '.*\\:',replacement = "",x = rownames(log2FC))
    GFF.Chr <- annot.GFF$V1
    GFF.Strand <- annot.GFF$V7
    GFF.Start <- annot.GFF$V4
    GFF.End <- annot.GFF$V5

    ### Now fill in log2FC data for later use in output
    # now assign Genes & Bins to log2FC
    log2FC$Genes <- annot.GFF$V10
    log2FC$Bins <- log2Bins

    # compute means and variances for the log2FC data frame
    log2FC$log2FC <- apply(log2(log2FC[,Indices2]+1),1,mean) - apply(log2(log2FC[,Indices1]+1),1,mean)
    log2FC$log2Variance <- apply(log2(log2FC[,Indices2]+1),1,var) + apply(log2(log2FC[,Indices1]+1),1,var)

    ### now we work on organizing the newlog2FC read count data
    # remove all exon bin rows with fewer than 2 mean reads across samples
    newlog2FC <- newlog2FC[which(apply(newlog2FC,1,mean) >= 2),]
    # now add EnsID to remove duplicated rows
    newlog2FC$EnsID <- sub(":.*","",rownames(newlog2FC))
    # identify duplicated rows
    dupRows <- which(duplicated(newlog2FC))
    # remove duplicated rows
    newlog2FC <- newlog2FC[-c(dupRows),]
    # remove EnsiD column
    newlog2FC <- newlog2FC[,-c(ncol(newlog2FC))]
    # set up newlog2FC "Bins"
    newBins <- sub(".*:","",rownames(newlog2FC))

    # take Ensembl ID numbers for easy parsing of newlog2FC data
    TempID <- rownames(newlog2FC)
    TempID <- gsub("\\:.*","",TempID)
    TempID <- gsub("\\..*","",TempID)
    TempID <- substr(TempID,(nchar(TempID[1]) - 9),nchar(TempID[1]))
    newIDs <- TempID
    rm(TempID)

    ### Add Ensembl ID column to annot.GFF data frame for later sorting
    TempID <- annot.GFF[,2]
    TempID <- gsub("\\:.*","",TempID)
    TempID <- gsub("\\..*","",TempID)
    TempID <- substr(TempID,(nchar(TempID[1]) - 9),nchar(TempID[1]))
    annot.GFF$EnsID <- TempID
    rm(TempID)

    ### remove ENSG00000 from ENSEMBL IDs to make them easier to sort (for looping through each unique ID)
    IDs <- gsub("\\:.*","",rownames(log2FC))
    IDs <- gsub("\\..*","",IDs)
    IDs <- substr(IDs,(nchar(IDs[1])-9),nchar(IDs[1]))

    ### Check to make that annot.GFF EnsID column and IDs vector are 100% identical
    EqualVectors <- all.equal(annot.GFF$EnsID,IDs)
    if (EqualVectors != TRUE){
        stop(call="The GFF file you provided does not match the GFF file used to count exonic reads. Please consider converting your GTF file to GFF again, and then re-counting exon reads.
    Otherwise, please make sure the same GFF file used to count exon reads is provided to ExCluster's input arguments.")
    }

    ###############################  IF TRUE, COMBINE COMMON EXONS  #####################################

    if (CombineExons == TRUE){
        ### ExCluster progression message
        cat("Combining common exons ...","",sep="\n")

        # Count number of rows in annot.GFF
        NROW_GFF <- nrow(annot.GFF)
        # Counter for counting new summed exon bins
        Counter <- 1
        # Start point for first newIDs gene row number
        newIDstart <- 1
        # Variable start will be used through loop
        Start <- 1
        # make newlog2FC list to dump results from CombineExons into
        temp_log2FC <- list()

        for (i in 2:nrow(annot.GFF)){
            ### check if we are on a new gene
            if ((IDs[i] != IDs[i-1]) == TRUE | i == NROW_GFF){
                # determine the Stop point of the gene
                if (i == NROW_GFF){
                    Stop <- i
                }else{
                    Stop <- (i-1)
                }
                # Rows of the log2FC data frame
                log2Rows <- Start:Stop

                ### now determine if the new gene is in newIDs
                if ((IDs[Start] == newIDs[newIDstart]) == TRUE){
                    # temporary newlog2FC data
                    TempIDs <- newIDs[newIDstart:(newIDstart+length(log2Rows))]
                    # determine the Stop point of newIDs
                    newIDstop <- newIDstart+length(which(TempIDs%in%IDs[Start]))-1
                    # Rows of the newIDs vector
                    newlog2Rows <- newIDstart:newIDstop
                    # grab all transcript combinations and take unique transcript combinations per gene
                    ExTrans <- annot.GFF[log2Rows,9]
                    # now grab only rows of log2FC matching & present in newlog2FC
                    matchingRows <- which(rownames(log2FC)[log2Rows]%in%rownames(newlog2FC)[newlog2Rows])
                    # now select only those elements from ExTrans based on the matchingRows
                    ExTrans <- ExTrans[c(matchingRows)]
                    # grab unique transcript combinations
                    uExTrans <- unique(ExTrans)
                    # grab original exon bins for the gene from log2FC
                    originalBins <- gsub(pattern = '.*\\:',replacement = "",x = rownames(log2FC)[log2Rows])
                    # now grab the new bins based on matching rows
                    newBins <- originalBins[matchingRows]

                    ### now run the combining exons script for this newlog2FC gene
                    if (length(uExTrans) > 1){
                        ## Make temporary read table for summing reads
                        TempReads <- newlog2FC[newlog2Rows,,drop=FALSE]
                        ### Make a matrix that will represent the presence or absence of transcripts
                        TransMatrix <- data.frame(matrix(data=0,nrow=length(uExTrans),ncol=ncol(TempReads)))
                        for (j in 1:length(uExTrans)){
                            uIndices <- which(ExTrans%in%uExTrans[j])
                            if (length(uIndices) > 1){
                                TransMatrix[j,] <- apply(TempReads[uIndices,],2,sum)
                                rownames(TransMatrix)[[j]] <- rownames(TempReads)[[uIndices[1]]]
                                # match new bins onto original log2FC rows
                                newBinIndices <- which(originalBins%in%newBins[uIndices])
                                # now make those original log2FC rows equal to the 1st newBins in uIndices
                                log2FC$Bins[newBinIndices] <- newBins[uIndices[1]]
                            }else{
                                TransMatrix[j,] <- TempReads[uIndices,]
                                rownames(TransMatrix)[[j]] <- rownames(TempReads)[[uIndices[1]]]
                                # match new bins onto original log2FC rows
                                newBinIndices <- which(originalBins%in%newBins[uIndices[1]])
                                # now make those original log2FC rows equal to the 1st newBins in uIndices
                                log2FC$Bins[newBinIndices] <- newBins[uIndices[1]]
                            }
                        }
                        # add the new TransMatrix to temp_log2FC if uExTrans > 1
                        temp_log2FC[[Counter]] <- TransMatrix
                        # increase the counter
                        Counter <- Counter + 1
                    }

                    ### now set up the newIDstart to match the next gene in the list, but not > length(newIDs)
                    newIDstart <- min((newIDstop+1),length(newIDs))
                }
                # now start the next gene
                Start <- i
            }
        }
        newlog2FC <- do.call(rbind,temp_log2FC)
    }
    # clean up annot.GFF
    rm(annot.GFF)
    # clean up
    rm(temp_log2FC)
    rm(log2Bins)
    rm(IDs)
    rm(newIDs)


    #####################################################################################################
    ###############################  COMPUTE MEANS, VARIANCES, ETC  #####################################
    #####################################################################################################

    ### ExCluster progression message
    cat("Pre-processing data ...","",sep="\n")

    ### now set up the log2Clusters data frame which will be used to run the analysis
    log2Clusters <- newlog2FC[,1:4]

    ### rename log2Clusters columns && set up additional columns
    colnames(log2Clusters) <- c("gene","type","log2FC","var")
    rownames(log2Clusters) <- rownames(log2Clusters)
    ### now that we have added the reads to log2Clusters, we can fill in more columns
    log2Clusters$gene <- gsub("\\:.*","",rownames(newlog2FC))
    log2Clusters$log2FC <- apply(log2(newlog2FC[,Indices2]+1),1,mean) - apply(log2(newlog2FC[,Indices1]+1),1,mean)
    log2Clusters$var <- apply(log2(newlog2FC[,Indices2]+1),1,var) + apply(log2(newlog2FC[,Indices1]+1),1,var)
    log2Clusters$meanreads <- apply(newlog2FC,1,mean)
    log2Clusters$clustnum <- NA
    log2Clusters$MaxClust <- 1
    log2Clusters$pval <- 1
    log2Clusters$padj <- 1

    ### now set add in read counts from newlog2FC
    # number of columns in log2Clusters
    NCols <- ncol(log2Clusters)
    # add in count data
    log2Clusters[,NCols:(NCols+ncol(newlog2FC)-1)] <- newlog2FC
    colnames(log2Clusters)[NCols:(NCols+ncol(newlog2FC)-1)] <- c(paste("cond1_rep",1:length(Indices1),sep=""),paste("cond2_rep",1:length(Indices2),sep=""))
    # Index which columns correspond to conditions 1 and 2 in the 'log2Clusters' data frame
    log2Indices1 <- grep("cond1_",colnames(log2Clusters))
    log2Indices2 <- grep("cond2_",colnames(log2Clusters))
    # now we can remove the 'newlog2FC' data frame
    rm(newlog2FC)

    ### Remove exons with 0 variance
    log2Clusters <- log2Clusters[which(log2Clusters$var > 0),]

    ### Remove exons with a log2 variance of > 16 (standard deviation of log2FC 4 or greater -- too high uncertainty)
    log2Clusters <- log2Clusters[which(log2Clusters$var <= 16),]

    ### First remove exon bins with fewer than 8 reads averaged between all replicates (cond1 & cond2)
    log2Clusters <- log2Clusters[which(log2Clusters$meanreads >= 5),]

    #### Remove genes with fewer than 4 reads on average in either condition
    GeneIDs <- unique(gsub("\\..*","",log2Clusters$gene))
    GeneIDs <- substr(GeneIDs,(nchar(GeneIDs)-9),nchar(GeneIDs))
    log2IDs <- gsub("\\..*","",log2Clusters$gene)
    log2IDs <- substr(log2IDs,(nchar(log2IDs)-9),nchar(log2IDs))

    Counter <- 1
    Start <- 1
    RemovalGenes <- NULL
    CounterG <- 1

    ### run loop to flag genes with fewer than 3 average reads in either condition 1 or condition 2
    for (i in 1:nrow(log2Clusters)){
        if ((GeneIDs[Counter] != log2IDs[i]) || i == nrow(log2Clusters)){
            ### determine where the stop is
            if (i == nrow(log2Clusters)){
                Stop <- i
            }else{
                Stop <- (i-1)
            }
            ### also remove genes with fewer than 4 average reads per gene in either condition
            Cond1_GeneExp <- mean(unlist(log2Clusters[Start:Stop,log2Indices1]))
            Cond2_GeneExp <- mean(unlist(log2Clusters[Start:Stop,log2Indices2]))
            Min_GeneExp <- min(Cond1_GeneExp,Cond2_GeneExp)
            if (Min_GeneExp < 2){
                RemovalGenes[CounterG] <- log2Clusters$gene[Start]
                CounterG <- CounterG + 1
            }
            Counter <- min((Counter+1),length(GeneIDs))
            Start <- i
        }
    }

    ### now grab indices for all genes flagged for removal (RemovalGenes) above\
    if (is.null(RemovalGenes) == FALSE){
        RemovalIndices <- which(log2Clusters$gene%in%RemovalGenes)
        ### remove these genes from log2Clusters
        log2Clusters <- log2Clusters[-c(RemovalIndices),]
    }

    ### Index log2Clusters row starst & stops for each gene -- greatly speeds up algorithm
    # Set up initial variables for indexing
    GeneIDs <- unique(gsub("\\..*","",log2Clusters$gene))
    GeneIDs <- substr(GeneIDs,(nchar(GeneIDs)-9),nchar(GeneIDs))
    log2IDs <- gsub("\\..*","",log2Clusters$gene)
    log2IDs <- substr(log2IDs,(nchar(log2IDs)-9),nchar(log2IDs))

    Start <- 1
    RowStart <- NULL
    RowEnd <- NULL
    Counter <- 1

    ### Run indexing loop
    for (i in 1:nrow(log2Clusters)){
        if ((GeneIDs[Counter] == log2IDs[i]) == FALSE || i == nrow(log2Clusters)){
            ## determine stop
            if (i == nrow(log2Clusters)){
                Stop <- i
            }else{
                Stop <- (i-1)
            }
            # now assign start/stop
            RowStart[Counter] <- Start
            RowEnd[Counter] <- Stop
            # now add 1 to counter & set "Start" to i
            Start <- i
            Counter <- Counter + 1
            if (i == nrow(log2Clusters) && Counter == length(GeneIDs)){
                # now assign start/stop
                RowStart[Counter] <- Start
                RowEnd[Counter] <- Stop
            }
        }
    }


    #####################################################################################################
    #############################  PREPARE & RUN CLUSTERING ALGORITHM  ##################################
    #####################################################################################################

    ### grab unique genes
    Genes <- unique(log2Clusters$gene)
    #### Count number of exons per gene
    log2Clusters$NumExons <- NA
    #### ultimately generate permutation tested pvals
    log2Clusters$pval <- NA
    ## do we run calcluations on this gene??
    ## to be used if nrows > 1 or if NbClust successfully runs
    RunClust <- FALSE
    ## random variables
    ClustCount <- 1
    NumClusters <- 1
    ## collect p-values per gene
    GenePVals <- array()
    ComputedPVals <- array()
    ## collect gene names for these arrays
    NamesGenePVals <- array()

    ### add columns with mean and variances for each condition in log2 space
    log2Clusters$mean1 <- apply(log2(log2Clusters[,log2Indices1]+1),1,mean)
    log2Clusters$mean2 <- apply(log2(log2Clusters[,log2Indices2]+1),1,mean)
    log2Clusters$var1 <- apply(log2(log2Clusters[,log2Indices1]+1),1,var)
    log2Clusters$var2 <- apply(log2(log2Clusters[,log2Indices2]+1),1,var)

    log2Clusters$var1 <- log2Clusters$var1 + 0.0000001
    log2Clusters$var2 <- log2Clusters$var2 + 0.0000001

    ######################################################################################################
    #########################################  NOW RUN EXCLUSTER  ########################################
    ######################################################################################################

    ### ExCluster progression message
    cat("Running main ExCluster module (very long) ...","",sep="\n")

    ## set up progress reporting on ExCluster
    ProgressPoints <- round(length(Genes)*c((1:10)/10))
    ProgressCounter <- 1

    ### loop through and hierarchically cluster all expressed genes with > 1 exon bin
    for (i in 1:length(Genes)){
        # display progress at every 10% of genes analyzed
        if (i == ProgressPoints[ProgressCounter]){
            cat(paste("ExCluster progress: ",(ProgressCounter*10),"%",sep=""),sep="\n")
            ProgressCounter <- min((ProgressCounter+1),10)
        }

        ## rows of log2Clusters datamatrix for Genes[i]
        rows <- RowStart[i]:RowEnd[i]
        ## number of rows in Genes[i]
        NRows <- length(rows)
        log2Clusters$NumExons[rows] <- NRows

        ### DM is the Data Matrix to make hellinger distance matrix
        Mean1 <- rep(log2Clusters$log2FC[rows],each=NRows)
        Var1 <- rep(log2Clusters$var[rows], each=NRows)
        Mean2 <- rep(log2Clusters$log2FC[rows],NRows)
        Var2 <- rep(log2Clusters$var[rows],NRows)

        ### run Effect Size computation
        ES_mat <- matrix(unlist(lapply(1:(NRows^2),FUN=Generate_ES_dists)),nrow=NRows,byrow=TRUE)
        Mean_ES <- mean(abs(unlist(ES_mat)))
        if (Mean_ES == 0){
            NRows <- 1
        }

        if (NRows > 3){
            HC <- hclust(dist(ES_mat,method = DISTANCE),method = LINKAGE)
            res <- Cutree_SD.Index()
            log2Clusters$clustnum[rows] <- res
        }else{
            log2Clusters$clustnum[rows] <- 1:length(rows)
        }
        ##Calc the number of clusters in the gene being tested
        NumClusters <- max(log2Clusters$clustnum[rows])
        log2Clusters$MaxClust[rows] <- max(log2Clusters$clustnum[rows])

        if (NumClusters > 1){
            ##Generate hclust object
            HC <- hclust(dist(ES_mat,DISTANCE),method=LINKAGE)
            ##cut tree
            HC_cutree <- cutree(HC,k=NumClusters)
            ## find row indices per cluster
            HC_indices <- lapply(1:NumClusters,FUN=Determine_HC_Distance)
            ## find the average euclidean distance between two most different clusters in real observed data
            Observed_clust_matrix <- matrix(unlist(lapply(1:NumClusters,FUN=Compare_HC_Distances)),nrow=NumClusters,ncol=NumClusters,byrow = TRUE)
            Observed_clust_dist <- max(apply(Observed_clust_matrix,1,sum)/(NumClusters-1))


            ###################  RUN PERMUATIONS PER GENE  ####################
            rm(HC_indices)
            rm(HC)
            rm(HC_cutree)

            ## grab overall log2FC of gene for simulation
            gene_log2FC <- log2(mean(apply(log2Clusters[rows,log2Indices2],2,sum))+1) - log2(mean(apply(log2Clusters[rows,log2Indices1],2,sum))+1)
            Gene.FC <- 2^gene_log2FC
            ## gene ratio 1
            Gene.R1 <- 1/(1+Gene.FC)
            ## gene ratio 2
            Gene.R2 <- Gene.FC/(1+Gene.FC)

            ### Run first 20 iterations
            Nullhypo_Dists <- generate_Nullhypo_Dists(NumIterations=20, gene_FC=Gene.FC, Gene_R1=Gene.R1, Gene_R2=Gene.R2, NumRows=NRows)
            # compute pval percentile in observed null hypothesis cluster distances
            f <- ecdf(Nullhypo_Dists)
            PVal <- 1 - f(Observed_clust_dist)

            ##############  NOW RUN 30 MORE ITERATIONS IF 'PVal' <= 0.5
            if (PVal <= 0.5){
                ### Run 30 more iterations
                Nullhypo_Dists_temp <- generate_Nullhypo_Dists(NumIterations=30, gene_FC=Gene.FC, Gene_R1=Gene.R1, Gene_R2=Gene.R2, NumRows=NRows)
                # combine these results with previous
                Nullhypo_Dists <- c(Nullhypo_Dists, Nullhypo_Dists_temp)
                # compute pval percentile in observed null hypothesis cluster distances
                f <- ecdf(Nullhypo_Dists)
                PVal <- 1 - f(Observed_clust_dist)
            }

            ##############  NOW RUN 50 MORE ITERATIONS IF 'PVal' <= 0.2
            if (PVal <= 0.2){
                ### Run 50 more iterations
                Nullhypo_Dists_temp <- generate_Nullhypo_Dists(NumIterations=50, gene_FC=Gene.FC, Gene_R1=Gene.R1, Gene_R2=Gene.R2, NumRows=NRows)
                # combine these results with previous
                Nullhypo_Dists <- c(Nullhypo_Dists, Nullhypo_Dists_temp)
                # compute pval percentile in observed null hypothesis cluster distances
                f <- ecdf(Nullhypo_Dists)
                PVal <- 1 - f(Observed_clust_dist)
            }

            ##############  NOW RUN 100 MORE ITERATIONS IF 'PVal' <= 0.05
            if (PVal <= 0.05){
                ### Run 100 more iterations
                Nullhypo_Dists_temp <- generate_Nullhypo_Dists(NumIterations=100, gene_FC=Gene.FC, Gene_R1=Gene.R1, Gene_R2=Gene.R2, NumRows=NRows)
                # combine these results with previous
                Nullhypo_Dists <- c(Nullhypo_Dists, Nullhypo_Dists_temp)
                # compute pval percentile in observed null hypothesis cluster distances
                f <- ecdf(Nullhypo_Dists)
                PVal <- 1 - f(Observed_clust_dist)
            }

            ##############  NOW RUN 300 MORE ITERATIONS IF 'PVal' <= 0.01
            if (PVal <= 0.01){
                ### Run 300 more iterations
                Nullhypo_Dists_temp <- generate_Nullhypo_Dists(NumIterations=300, gene_FC=Gene.FC, Gene_R1=Gene.R1, Gene_R2=Gene.R2, NumRows=NRows)
                # combine these results with previous
                Nullhypo_Dists <- c(Nullhypo_Dists, Nullhypo_Dists_temp)
                # compute pval percentile in observed null hypothesis cluster distances
                f <- ecdf(Nullhypo_Dists)
                PVal <- 1 - f(Observed_clust_dist)
            }

            ##############  NOW RUN 500 MORE ITERATIONS IF 'PVal' <= 0.004
            if (PVal <= 0.004){
                ### Run 500 more iterations
                Nullhypo_Dists_temp <- generate_Nullhypo_Dists(NumIterations=500, gene_FC=Gene.FC, Gene_R1=Gene.R1, Gene_R2=Gene.R2, NumRows=NRows)
                # combine these results with previous
                Nullhypo_Dists <- c(Nullhypo_Dists, Nullhypo_Dists_temp)
                # compute pval percentile in observed null hypothesis cluster distances
                f <- ecdf(Nullhypo_Dists)
                PVal <- 1 - f(Observed_clust_dist)
            }
            ### Now if PVal == 0, use gamma distributions to estimate PValue better, but it can't be higher than 0.0001
            if (PVal == 0){
                # p-value based on gamma function
                PVal <- gm_mean(pgamma(Observed_clust_dist, Nullhypo_Dists, lower = FALSE) * gamma(Nullhypo_Dists))
                PVal <- min(PVal, 0.00001)
            }

            # clean up
            rm(Nullhypo_Dists)
            rm(Observed_clust_dist)
            # assign pvals & Gene Names for stats
            log2Clusters$pval[rows] <- PVal
            GenePVals[i] <- PVal
            NamesGenePVals[i] <- log2Clusters$gene[rows[1]]
        }else{
            GenePVals[i] <- NA
            NamesGenePVals[i] <- log2Clusters$gene[rows[1]]
        }
        rm(NRows)
    }

    ######################################## STATISTICS ########################################

    ### ExCluster progression message
    cat("Running statistical tests ...","",sep="\n")

    ### data frame for testing p-values
    Catch_PVals <- data.frame(NamesGenePVals,GenePVals)
    colnames(Catch_PVals) <- c("EnsG","pvals")
    # remove NAs
    Catch_PVals <- Catch_PVals[which(complete.cases(Catch_PVals[,2])),]
    # sort PVals from lowest to highest
    Catch_PVals <- Catch_PVals[order(Catch_PVals[,2]),]

    if (nrow(Catch_PVals) >= 4000){
        ### normalize p-values
        # estimate the number of 'truly spliced genes'
        pvalCutoff <- EstNullHypo(Catch_PVals[,2])
        # estimate the number of truly spliced genes
        NumTrue <- nrow(Catch_PVals[Catch_PVals[,2] <= pvalCutoff,])
        # estimate the number of null hypothesis genes
        NumNull <- nrow(Catch_PVals[Catch_PVals[,2] > pvalCutoff,])
        ### set the minimum value a null hypothesis p-value can be corrected to
        minNullPVal <- min(runif(NumNull/(pvalCutoff/50),0,1))

        ### if we have > 0 NumTrue p-values
        if (NumTrue > 0){
            # Now give the null hypothesis p-values uniform distriubtions
            Catch_PVals[(NumTrue+1):(NumTrue+NumNull),2] <- sort(runif(n=NumNull,minNullPVal,1))
            # Now divide all true p-values by a factor based on the lowest null hypothesis value
            Catch_PVals[1:(NumTrue),2] <- Catch_PVals[c(1:NumTrue),2]/(Catch_PVals[NumTrue,2]/minNullPVal)
            # Now compute FDR values
            FDRs <- FDRcalc(Catch_PVals[,2],NumTrue,NumNull)
        }else{

            ### if we have no truly significant p-values from our estimate
            Catch_PVals[,2] <- sort(runif(n=NumNull,0,1))
            FDRs <- rep(1,nrow(Catch_PVals))
        }

    }else{

        ### This part of the code will be run if less than 4000 expressed genes were identified
        warning(call = "You have entered transcriptome data which contains fewer than 4000 expressed genes. ExCluster cannot run advanced statistics (yielding better tuned FDRs) with < 4000 expressed genes.
Please avoid running ExCluster on sub-sections of a genome (i.e. only certain chromosomes), as background expressed genes are used to tune statistics.

If you are running test data to get ExCluster working, ignore this warning.")

        # because we have too few pvalues (< 4000), we cannot run our advanced stats
        # we therefore run a basic correction with Benjamini-Hochberg
        pvalCutoff <- 0.05
        # estimate the number of truly spliced genes
        NumTrue <- nrow(Catch_PVals[Catch_PVals[,2] <= pvalCutoff,])
        # estimate the number of null hypothesis genes
        NumNull <- nrow(Catch_PVals[Catch_PVals[,2] > pvalCutoff,])
        # set the minimum value a null hypothesis p-value can be corrected to
        minNullPVal <- min(runif(NumNull/2,0,1))

        ### if we have truly differentially spliced genes
        if (NumTrue > 0){
            # Now give the null hypothesis p-values uniform distriubtions
            Catch_PVals[(NumTrue+1):(NumTrue+NumNull),2] <- sort(runif(n=NumNull,minNullPVal,1))
            # Now divide all true p-values by a factor based on the lowest null hypothesis value
            Catch_PVals[1:(NumTrue),2] <- Catch_PVals[c(1:NumTrue),2]/(Catch_PVals[NumTrue,2]/minNullPVal)
            # Now compute FDR values
            FDRs <- p.adjust(p = Catch_PVals[,2],method="BH")
        }else{

            ### if we have no truly significant p-values from our estimate
            Catch_PVals[,2] <- sort(runif(n=NumNull,0,1))
            FDRs <- rep(1,nrow(Catch_PVals))
        }
    }

    # Make statistics table per gene with pvals and FDRs
    Stats_table <- data.frame(Catch_PVals,FDRs)

    ######################################################################################################
    ###############################  FINAL DATA PREPARATION FOR OUTPUT ###################################
    ######################################################################################################

    ### ExCluster progression message
    cat("Final ExCluster data organization & output ...","",sep="\n")

    # set up statistics & cluster columns in log2FC
    log2FC$pval <- NA
    log2FC$FDR <- NA
    log2FC$Cluster <- NA

    ### remove ENSG00000 from ENSEMBL IDs to make them easier to sort (for looping through each unique ID)
    # do this to add an IDs column to log2FC
    IDs <- gsub("\\:.*","",rownames(log2FC))
    IDs <- gsub("\\..*","",IDs)
    IDs <- substr(IDs,(nchar(IDs[1])-9),nchar(IDs[1]))
    # Now add an IDs column to the log2FC data frame
    log2FC$IDs <- IDs

    # now do the same for the IDs in the stats table
    IDs <- gsub("\\:.*","",Stats_table[,1])
    IDs <- gsub("\\..*","",IDs)
    IDs <- substr(IDs,(nchar(IDs[1])-9),nchar(IDs[1]))
    # Now add the IDs to Stats_table
    Stats_table$IDs <- IDs

    # now do the same for IDs in the log2Clusters
    IDs <- gsub("\\:.*","",log2Clusters$gene)
    IDs <- gsub("\\..*","",IDs)
    IDs <- substr(IDs,(nchar(IDs[1])-9),nchar(IDs[1]))
    # Now add the IDs to log2Clusters
    log2Clusters$IDs <- IDs
    # Also add the exon bin number to log2Clusters for matching
    log2Clusters$Bin <- gsub(pattern = '.*\\:',replacement = "",x = rownames(log2Clusters))

    # now order on EnsG
    log2FC <- log2FC[order(rownames(log2FC)),]
    Stats_table <- Stats_table[order(Stats_table$EnsG),]
    log2Clusters <- log2Clusters[order(rownames(log2Clusters)),]

    # set up counter for Stats_table genes
    StatsCounter <- 1
    # set up a Start location counter for the log2Clusters dataframe
    log2Start <- 1
    # Start location for first gene
    Start <- 1

    ### now run through log2Clusters & Stats_table and match gene data, transfering data to log2F for output
    for (i in 2:nrow(log2FC)){

        ### First determine if we have ended the current gene of interest
        if ((log2FC$IDs[i] != log2FC$IDs[(i-1)]) == TRUE | i == nrow(log2FC)){
            # Now determine the stop location for the current gene
            if (i == nrow(log2FC)){
                Stop <- i
            }else{
                Stop <- (i-1)
            }

            ### Now check if the Stats_table ID is equal to the current log2FC gene ID
            if (log2FC$IDs[Start] == Stats_table$IDs[StatsCounter]){
                # now assign stats data to the log2FC data frame, on rows Start:Stop
                log2FC$pval[Start:Stop] <- Stats_table$pvals[StatsCounter]
                log2FC$FDR[Start:Stop] <- Stats_table$FDRs[StatsCounter]
                StatsCounter <- min((StatsCounter+1), nrow(Stats_table))
            }

            ### Now check if the log2Clusters ID is equal to the current log2FC gene ID
            if (log2FC$IDs[Start] == log2Clusters$IDs[log2Start]){
                # if above is TRUE, grab data for the next 300 rows in log2Clusters
                # make sure we don't run out of rows in log2Start
                Upper <- min(300,(nrow(log2Clusters)-log2Start))
                TempData <- log2Clusters[log2Start:(log2Start+Upper),]
                # isolate data for only the current gene
                TempData <- TempData[which(TempData$IDs%in%TempData$IDs[1]),]
                # grab the exon bins from this gene
                ExonBins <- TempData$Bin
                # now loop through the exon bins in TempData & assign cluster numbers to log2FC
                for (j in 1:length(ExonBins)){
                    BinIndices <- which(log2FC$Bins[Start:Stop]%in%ExonBins[j])
                    log2FC$Cluster[Start:Stop][BinIndices] <- TempData$clustnum[j]
                }
                # Now set log2Start to the beginning of the next log2Clusters gene
                log2Start <- min(log2Start+length(ExonBins), nrow(log2Clusters))
            }

            ### Reset the gene start for the next log2FC gene
            Start <- i
        }
    }
    # clean up
    rm(log2Clusters)
    rm(Stats_table)

    # separate EnsG and exon bin number
    ExonBins <- gsub(pattern = '.*\\:',replacement = "",x = rownames(log2FC))
    EnsG <- gsub(pattern = '\\:.*',replacement = "",x = rownames(log2FC))

    ### now re-arrange log2FC columns to be more logical
    final.log2FC <- data.frame(EnsG, ExonBins ,log2FC$Genes,log2FC$Cluster,log2FC$log2FC,log2FC$log2Variance,log2FC$pval,log2FC$FDR,GFF.Chr,GFF.Strand,GFF.Start,GFF.End)
    # rename columns
    colnames(final.log2FC) <- c("EnsG","Exon_bin","Gene_name","Cluster","log2FC","log2Var","pval","FDR","Chr","Strand","Start","End")
    # re-add the read counts at the end of the dataframe
    final.log2FC <- data.frame(final.log2FC,log2FC[,1:length(cond.Nums)])
    # clean up log2FC
    rm(log2FC)

    ### if out.Dir was specified, write out ExCluster results
    if (is.null(out.Dir) == FALSE){
        # check to make sure out.Dir path can be written to
        WriteCheck <- file.access(out.Dir, mode=2)
        if (WriteCheck == 0){
            # now check to see if a result.Filename was specified
            if (is.null(result.Filename) == TRUE){
                # assign full filename
                fullFilename <- paste(out.Dir,"/ExClust_Results",sep="")
            }else{
                fullFilename <- paste(out.Dir,"/",result.Filename,sep="")
            }
            # make sure we're not over-writing an existing file
            if (file.exists(paste(fullFilename,".txt",sep="")) == TRUE){
                fullFilename <- paste(fullFilename,Sys.time(),sep="")
            }
            write.table(final.log2FC,file=paste(fullFilename,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
        }else{
            print("Error! You specified an out.Dir directory for ExCluster that could not be written to. Please provide a valid directory that ExCluster can write files to. ExCluster will finish running, however your file has not been written.")
        }
    }

    ### check to see if we plot the results from ExCluster
    if (plot.Results == TRUE){
        if (is.null(out.Dir) == FALSE){
            # new images folder to write to
            out.Dir <- paste(out.Dir,"/exon_log2FC_images/",sep="")
            # create the images directly if it doesn't already exist
            WriteCheck <- file.access(out.Dir, mode=2)
            if (WriteCheck != 0){
                dir.create(path = out.Dir,recursive = TRUE,)
            }
            # make sure that FDR.cutoff is less than 0.2, or set it back to 0.05
            if (FDR.cutoff > 0.2){
                warning(call="FDR.cutoff was assigned higher than 0.2, which is not allowed. Using FDR cutoffs above 20% generates unnecessarily high false discovery rates. Please set FDR.cutoff to less than or equal to 0.2 in the future.")
                FDR.cutoff <- 0.05
            }
            # run plotting function
            plotExClustResults(results.Data=final.log2FC, out.Dir=out.Dir, FDR.cutoff=FDR.cutoff)
        }else{
            print("Error! You specified that plot.Results=TRUE, but you did not specify an directory for output. out.Dir must be assigned if plot.Results = TRUE -- your exon log2FCs have not been plotted.")
            print("However, if you have saved your ExCluster results to a variable or file, you may run the plotExClustResults function directly -- please consult the vignette.")
        }
    }
    ### return final data frame resutls
    return(final.log2FC)
}
