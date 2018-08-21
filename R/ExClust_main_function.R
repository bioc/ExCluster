ExClust_main_function <- function(newlog2FC=NULL, Indices1=NULL, Indices2=NULL, NumReps1=NULL, NumReps2=NULL){

    ################################  SETUP PARAMETERS & FUNCTIONS ################################

    ## clustering method = euclidean
    DISTANCE <- "euclidean"

    ## linkage == complete
    LINKAGE <- "complete"

    ### Function for running null hypothesis non-parametric iterations
    generate_Nullhypo_Dists <- function(NumIterations=NULL,gene_FC=NULL,Gene_R1=NULL,Gene_R2=NULL, NumRows=NULL){
        ## rows = genes, columns = exons
        Sim_log2FC <- matrix(NA,nrow=NumIterations,ncol=NumRows)
        Sim_log2var <- matrix(NA,nrow=NumIterations,ncol=NumRows)

        for (l in seq(NumRows)){
            # generate rnorm vector of random null hypothesis sampling
            SimulatedMatrix1 <- abs(rnorm(n = NumIterations*(length(log2Indices1)), sd=sqrt(log2Clusters$var1[rows[l]]),
                                          mean=log2(((apply(log2Clusters[rows[l],c(log2Indices1,log2Indices2)],1,mean)*2)*Gene_R1)+1)))
            SimulatedMatrix2 <- abs(rnorm(n = NumIterations*(length(log2Indices2)), sd=sqrt(log2Clusters$var2[rows[l]]),
                                          mean=log2(((apply(log2Clusters[rows[l],c(log2Indices1,log2Indices2)],1,mean)*2)*Gene_R2)+1)))
            # transform these vectors into matrices
            SimulatedMatrix1 <- matrix(SimulatedMatrix1, nrow = NumIterations)
            SimulatedMatrix2 <- matrix(SimulatedMatrix2, nrow = NumIterations)

            # compute log2FC means and variances from above matrices
            Sim_log2FC[,l] <- apply(SimulatedMatrix2,1,mean) - apply(SimulatedMatrix1,1,mean)
            Sim_log2var[,l] <- apply(SimulatedMatrix2,1,var) + apply(SimulatedMatrix1,1,var)
        }

        Nonparam_Clust_Dists <-vapply(seq(NumIterations),Permuted_ES_nullhypo,Sim_log2FC=Sim_log2FC,
                                      DISTANCE=DISTANCE, LINKAGE=LINKAGE, Sim_log2var=Sim_log2var,
                                      NumClusters=NumClusters, NumRows=NumRows,FUN.VALUE = c(Res=0))
        return(Nonparam_Clust_Dists)
    }

    ### ExCluster progression message
    cat("Pre-processing data ...","",sep="\n")

    ### now set up the log2Clusters data frame which will be used to run the analysis
    log2Clusters <- newlog2FC[,seq(4)]

    ### rename log2Clusters columns & rows
    colnames(log2Clusters) <- c("gene","type","log2FC","var")
    rownames(log2Clusters) <- rownames(log2Clusters)
    ### now that we have added the reads to log2Clusters, we can fill in more columns
    # gene_id values without exon bins
    log2Clusters$gene <- gsub("\\:.*","",rownames(newlog2FC))
    # log2FC values
    log2Clusters$log2FC <- apply(log2(newlog2FC[,Indices2]+1),1,mean) -
        apply(log2(newlog2FC[,Indices1]+1),1,mean)
    # log2 variances calculations per condition
    log2Vars1 <- apply(log2(newlog2FC[,Indices1]+1),1,var)
    log2Vars2 <- apply(log2(newlog2FC[,Indices2]+1),1,var)
    # the log2FC variance for log2Clusters is the sum of the previous 2 variances
    log2Clusters$var <- log2Vars1 + log2Vars2
    # mean reads across conditions (normal space, not log2 space)
    log2Clusters$meanreads <- apply(newlog2FC,1,mean)
    # additional metadata columns to be used later
    log2Clusters$clustnum <- NA
    log2Clusters$MaxClust <- 1
    log2Clusters$pval <- 1
    log2Clusters$padj <- 1

    ### now set add in read counts from newlog2FC
    # number of columns in log2Clusters
    NCols <- ncol(log2Clusters)
    # add in count data
    log2Clusters[,seq(NCols,(NCols+ncol(newlog2FC)-1))] <- newlog2FC
    colnames(log2Clusters)[seq(NCols,(NCols+ncol(newlog2FC)-1))] <-
        c(paste("cond1_rep",seq(length(Indices1)),sep=""),paste("cond2_rep",seq(length(Indices2)),sep=""))
    # Index which columns correspond to conditions 1 and 2 in the 'log2Clusters' data frame
    log2Indices1 <- grep("cond1_",colnames(log2Clusters))
    log2Indices2 <- grep("cond2_",colnames(log2Clusters))
    # now we can remove the 'newlog2FC' data frame
    rm(newlog2FC,envir=sys.frame(-1))

    ### now we remove exons with a log2 read variance of > 9 in either condition (very high variances)
    log2Clusters <- log2Clusters[which(log2Vars1 <= 9 & log2Vars2 <= 9),]
    # clean up
    rm(log2Vars1)
    rm(log2Vars2)

    ### Remove exons with 0 variance
    log2Clusters <- log2Clusters[which(log2Clusters$var > 0),]

    ### First remove exon bins with fewer than 4 max reads across all conditions/replicates, for a given exon
    log2Clusters <- log2Clusters[which(apply(log2Clusters[,c(log2Indices1,log2Indices2)],1,max) >= 4),]

    #### Remove genes with fewer than 4 reads on average in either condition
    log2IDs <- gsub("\\..*","",log2Clusters$gene)
    log2IDs <- substr(log2IDs,(nchar(log2IDs)-9),nchar(log2IDs))
    Counter <- 1
    Start <- 1
    RemovalIndices <- vector("list", 1)

    ### run loop to flag genes with fewer than 4 reads in all exons of condition 1 or 2
    for (i in seq(nrow(log2Clusters))){
        nextValue <- min((i+1),nrow(log2Clusters))
        if ((log2IDs[i] != log2IDs[nextValue]) || i == nrow(log2Clusters)){
            # stop at the current gene
            Stop <- i
            # compute the maximum reads in each condition
            Cond1_GeneExp <- max(apply(log2Clusters[seq(Start,Stop),log2Indices1],1,max))
            Cond2_GeneExp <- max(apply(log2Clusters[seq(Start,Stop),log2Indices2],1,max))
            # take the minimum of these 2
            Min_GeneExp <- min(Cond1_GeneExp,Cond2_GeneExp)
            # if the minimum is less than 8 reads, remove this gene
            if (Min_GeneExp < 4){
                RemovalIndices[[Counter]] <- seq(Start,Stop)
                Counter <- Counter + 1
            }
            Start <- i+1
        }
    }

    ### now grab indices for all genes flagged for removal (RemovalIndices) above\
    if (length(unlist(RemovalIndices)) > 0){
        RemovalIndices <- unlist(RemovalIndices)
        ### remove these genes from log2Clusters
        log2Clusters <- log2Clusters[-c(RemovalIndices),]
    }

    ### Index log2Clusters row starst & stops for each gene -- greatly speeds up algorithm
    # Set up initial variables for indexing
    log2IDs <- gsub("\\..*","",log2Clusters$gene)
    log2IDs <- substr(log2IDs,(nchar(log2IDs)-9),nchar(log2IDs))

    Start <- 1
    RowStart <- NULL
    RowEnd <- NULL
    Counter <- 1

    ### Run indexing loop
    for (i in seq(nrow(log2Clusters))){
        nextValue <- min((i+1),nrow(log2Clusters))
        if ((log2IDs[i] != log2IDs[nextValue]) || i == nrow(log2Clusters)){
            # stop = ith row
            Stop <- i
            # now assign start/stop
            RowStart[Counter] <- Start
            RowEnd[Counter] <- Stop
            # now add 1 to counter & set "Start" to i
            Start <- (i+1)
            Counter <- Counter + 1
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

    ### now make any row with log2 mean < 0.33 & log2 variance < 0.33 == c(1,0,0)
    # the purpose of this is to eliminate zero read counts & zero variances that may persist
    # do this for condition 1 first
    # grab the total number of cells in the matrix
    num.Cells1 <- length(unlist(log2Clusters[which(log2Clusters$mean1 <= 0.33 & log2Clusters$var1 <= 0.33), log2Indices1]))
    # now make these cells either c(1,0,1) or c(1,0,0)
    if (num.Cells1 > NumReps1){
        log2Clusters[which(log2Clusters$mean1 <= 0.33 & log2Clusters$var1 <= 0.33), log2Indices1] <-
            c(rep(1,num.Cells1/NumReps1),rep(0,num.Cells1/NumReps1),sample(c(0,1),size=(num.Cells1-(2*(num.Cells1/NumReps1))),replace=TRUE))
    }
    # now do the same for condition 2
    num.Cells2 <- length(unlist(log2Clusters[which(log2Clusters$mean2 <= 0.33 & log2Clusters$var2 <= 0.33), log2Indices2]))
    if (num.Cells2 > NumReps2){
        log2Clusters[which(log2Clusters$mean2 <= 0.33 & log2Clusters$var2 <= 0.33), log2Indices2] <-
            c(rep(1,num.Cells2/NumReps2),rep(0,num.Cells2/NumReps2),sample(c(0,1),size=(num.Cells2-(2*(num.Cells2/NumReps2))),replace=TRUE))
    }

    ### now we re-compute log2 means and variances for each condition
    log2Clusters$mean1 <- apply(log2(log2Clusters[,log2Indices1]+1),1,mean)
    log2Clusters$mean2 <- apply(log2(log2Clusters[,log2Indices2]+1),1,mean)
    log2Clusters$var1 <- apply(log2(log2Clusters[,log2Indices1]+1),1,var)
    log2Clusters$var2 <- apply(log2(log2Clusters[,log2Indices2]+1),1,var)

    ### now make sure all genes with < 2 mean log2 read have minimum variance of 0.33
    if (length(log2Clusters$var1[which(log2Clusters$mean1 < 2)]) > 0){
        log2Clusters$var1[which(log2Clusters$mean1 < 2)] <- 0.33
    }
    if (length(log2Clusters$var2[which(log2Clusters$mean2 < 2)]) > 0){
        log2Clusters$var2[which(log2Clusters$mean2 < 2)] <- 0.33
    }

    ### now we make sure all log2 variances are at least 1e-08 (very low, no variance should be this low)
    log2Clusters$var1 <- log2Clusters$var1 + 0.00000001
    log2Clusters$var2 <- log2Clusters$var2 + 0.00000001

    ### lastly we re-compute the log2FC and log2var for the gene
    log2Clusters$log2FC <- log2Clusters$mean2 - log2Clusters$mean1
    log2Clusters$var <- log2Clusters$var1 + log2Clusters$var2

    ######################################################################################################
    #########################################  NOW RUN EXCLUSTER  ########################################
    ######################################################################################################

    ### ExCluster progression message
    cat("Running main ExCluster module (very long) ...","",sep="\n")

    ## set up progress reporting on ExCluster
    ProgressPoints <- round(length(Genes)*c((seq(10))/10))
    ProgressCounter <- 1

    ### loop through and hierarchically cluster all expressed genes with > 1 exon bin
    for (i in seq(length(Genes))){
        # display progress at every 10% of genes analyzed
        if (i == ProgressPoints[ProgressCounter]){
            cat(paste("ExCluster progress: ",(ProgressCounter*10),"%",sep=""),sep="\n")
            ProgressCounter <- min((ProgressCounter+1),10)
        }

        # rows of log2Clusters datamatrix for Genes[i]
        rows <- seq(RowStart[i],RowEnd[i])
        # number of rows in Genes[i]
        NRows <- length(rows)
        # number of exons in this gene
        log2Clusters$NumExons[rows] <- NRows

        if (NRows > 1){
            ### Now grab the log2FC means and variances for each exon to make Effect Size difference matrices
            Mean1 <- as.numeric(log2Clusters$log2FC[rows])
            Var1 <- as.numeric(log2Clusters$var[rows])

            ### run Effect Size computation
            ES_mat <- Generate_ES_dists(Mean1,Var1)
            Mean_ES <- mean(abs(unlist(ES_mat)))
            if (Mean_ES == 0){
                NRows <- 1
            }

            if (NRows > 3){
                HC <- hclust(dist(ES_mat,method = DISTANCE),method = LINKAGE)
                HC_cutree <- Cutree_SD.Index(DISTANCE, NRows, HC, ES_mat)
                log2Clusters$clustnum[rows] <- HC_cutree
            }else{
                log2Clusters$clustnum[rows] <- seq(length(rows))
            }
        }else{
            log2Clusters$clustnum[rows] <- 1
        }

        ##Calc the number of clusters in the gene being tested
        NumClusters <- max(log2Clusters$clustnum[rows])
        log2Clusters$MaxClust[rows] <- max(log2Clusters$clustnum[rows])

        if (NumClusters > 1 && NRows > 1){
            ##Generate hclust object
            ES_Dists <- as.matrix(dist(ES_mat,method = DISTANCE))
            HC <- hclust(as.dist(ES_Dists),method=LINKAGE)
            ##cut tree
            HC_cutree <- cutree(HC,k=NumClusters)
            ## find row indices per cluster
            HC_indices <- lapply(seq(NumClusters),FUN=Determine_HC_Distance, HC_cutree=HC_cutree)
            ## find the average euclidean distance between two most different clusters in real observed data
            Observed_clust_dist <- Compare_HC_Distances(NumClusters, ES_Dists, HC_indices)
            # clean up
            rm(HC_indices)
            rm(HC)
            rm(HC_cutree)
            rm(ES_Dists)
            rm(ES_mat)

            ###################  RUN PERMUATIONS PER GENE  ####################

            ## grab overall log2FC of gene for simulation
            gene_log2FC <- log2(mean(apply(log2Clusters[rows,log2Indices2],2,sum))+1) -
                log2(mean(apply(log2Clusters[rows,log2Indices1],2,sum))+1)
            Gene.FC <- 2^gene_log2FC
            ## gene ratio 1
            Gene.R1 <- 1/(1+Gene.FC)
            ## gene ratio 2
            Gene.R2 <- Gene.FC/(1+Gene.FC)

            ### Run first 20 iterations
            Nullhypo_Dists <- generate_Nullhypo_Dists(NumIterations=20, gene_FC=Gene.FC,
                                                      Gene_R1=Gene.R1, Gene_R2=Gene.R2, NumRows=NRows)
            # compute pval percentile in observed null hypothesis cluster distances
            f <- ecdf(Nullhypo_Dists)
            PVal <- 1 - f(Observed_clust_dist)

            ##############  NOW RUN 30 MORE ITERATIONS IF 'PVal' <= 0.5
            if (PVal < 0.5){
                ### Run 30 more iterations
                Nullhypo_Dists_temp <- generate_Nullhypo_Dists(NumIterations=30, gene_FC=Gene.FC,
                                                               Gene_R1=Gene.R1, Gene_R2=Gene.R2, NumRows=NRows)
                # combine these results with previous
                Nullhypo_Dists <- c(Nullhypo_Dists, Nullhypo_Dists_temp)
                # compute pval percentile in observed null hypothesis cluster distances
                f <- ecdf(Nullhypo_Dists)
                PVal <- 1 - f(Observed_clust_dist)
            }

            ##############  NOW RUN 50 MORE ITERATIONS IF 'PVal' <= 0.2
            if (PVal < 0.2){
                ### Run 50 more iterations
                Nullhypo_Dists_temp <- generate_Nullhypo_Dists(NumIterations=50, gene_FC=Gene.FC,
                                                               Gene_R1=Gene.R1, Gene_R2=Gene.R2, NumRows=NRows)
                # combine these results with previous
                Nullhypo_Dists <- c(Nullhypo_Dists, Nullhypo_Dists_temp)
                # compute pval percentile in observed null hypothesis cluster distances
                f <- ecdf(Nullhypo_Dists)
                PVal <- 1 - f(Observed_clust_dist)
            }

            ##############  NOW RUN 100 MORE ITERATIONS IF 'PVal' <= 0.05
            if (PVal < 0.05){
                ### Run 100 more iterations
                Nullhypo_Dists_temp <- generate_Nullhypo_Dists(NumIterations=100, gene_FC=Gene.FC,
                                                               Gene_R1=Gene.R1, Gene_R2=Gene.R2, NumRows=NRows)
                # combine these results with previous
                Nullhypo_Dists <- c(Nullhypo_Dists, Nullhypo_Dists_temp)
                # compute pval percentile in observed null hypothesis cluster distances
                f <- ecdf(Nullhypo_Dists)
                PVal <- 1 - f(Observed_clust_dist)
            }

            ##############  NOW RUN 300 MORE ITERATIONS IF 'PVal' <= 0.01
            if (PVal < 0.01){
                ### Run 300 more iterations
                Nullhypo_Dists_temp <- generate_Nullhypo_Dists(NumIterations=300, gene_FC=Gene.FC,
                                                               Gene_R1=Gene.R1, Gene_R2=Gene.R2, NumRows=NRows)
                # combine these results with previous
                Nullhypo_Dists <- c(Nullhypo_Dists, Nullhypo_Dists_temp)
                # compute pval percentile in observed null hypothesis cluster distances
                f <- ecdf(Nullhypo_Dists)
                PVal <- 1 - f(Observed_clust_dist)
            }

            ##############  NOW RUN 500 MORE ITERATIONS IF 'PVal' <= 0.004
            if (PVal < 0.004){
                ### Run 500 more iterations
                Nullhypo_Dists_temp <- generate_Nullhypo_Dists(NumIterations=500, gene_FC=Gene.FC,
                                                               Gene_R1=Gene.R1, Gene_R2=Gene.R2, NumRows=NRows)
                # combine these results with previous
                Nullhypo_Dists <- c(Nullhypo_Dists, Nullhypo_Dists_temp)
                # compute pval percentile in observed null hypothesis cluster distances
                f <- ecdf(Nullhypo_Dists)
                PVal <- 1 - f(Observed_clust_dist)
            }
            ### Now if PVal == 0, use gamma distributions to estimate PValue better, but it can't be higher than 0.0001
            if (PVal == 0){
                # p-value based on gamma function
                PVal <- gm_mean(pgamma(Observed_clust_dist, Nullhypo_Dists, lower.tail = FALSE) * gamma(Nullhypo_Dists))
                PVal <- min(PVal, 0.0001)
            }

            # clean up
            rm(Nullhypo_Dists)
            rm(Observed_clust_dist)
            # assign pvals & Gene Names for stats
            log2Clusters$pval[rows] <- PVal
        }
        rm(NRows)
    }
    return(log2Clusters)
}
