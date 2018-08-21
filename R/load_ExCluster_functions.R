############################################################################
################################ GFF_convert ###############################
############################################################################

### take rtracklayer objects & reformat them into standard GTF
reformat_GTF <- function(rtracklayer.GTF=NULL){

    ### This function converts rtracklayer imported GTF files into flattened GFF files
    ### it functions by converting rtracklayer format back to GTF, and then calls GFF_convert()

    ### export rtracklayer object to GTF file (temporary gzip file)
    export(rtracklayer.GTF, format = "GTF", con=paste(tempdir(),"/temp.gtf",sep=""))

    ### now read this file back in
    annot.GTF <- read.table(file=paste(tempdir(),"/temp.gtf",sep=""),sep="\t", stringsAsFactors=FALSE)

    ### clean up
    file.remove(paste(tempdir(),"/temp.gtf",sep=""))

    ### return result & end function
    return(annot.GTF)
}

# This function takes the input of GTF data in the specific format produced in the GFF_collapse() function
# Normally this would be GTF data for a single gene
CollapseExons <- function(GR.data){
    # number of rows in the GTF dataframe to be made
    NRows <- nrow(GR.data)
    if (NRows > 1){
        # generate GRanges for GTF
        GR.gtf <- GRanges(seqnames=GR.data$V1, ranges=IRanges(GR.data$V4,GR.data$V5),tx_id=GR.data$tx.id)
        # collapse GRanges for non-overlapping exon bins for GFF (loses metadata)
        GR.gff <- disjoin(GR.gtf)
        # assuming NRows was > 1, now recount the new NRows in the GFF file
        NRows <- length(GR.gff)
        # now map indices for which rows of GR.gff match which rows in GR.gtf
        indices.gtf <- as(findOverlaps(GR.gff, GR.gtf), "List")
        GFF.transcripts <- NULL
        for (n in seq(NRows)){
            GFF.transcripts[n] <- paste(GR.gtf$tx_id[indices.gtf[[n]]],collapse='+')
        }
        GFF.start <- start(GR.gff)
        GFF.end <- end(GR.gff)
    }else{
        GFF.transcripts <- GR.data$tx.id[1]
        GFF.start <- GR.data$V4[1]
        GFF.end <- GR.data$V5[1]
    }

    # now make the rest of the columns for data frame output of GFF data
    GFF.chr <- rep(GR.data$V1[1],NRows)
    GFF.geneid <- rep(GR.data$gene.id[1],NRows)
    GFF.names <- rep(GR.data$gene.name[1],NRows)

    ### generate the data frame
    GFF.dataframe <- cbind(GFF.chr,GFF.geneid,GFF.start,GFF.end,GFF.transcripts,GFF.names)
    # return the data
    return(GFF.dataframe)
}

### this function indexes the start/stop rows of each gene in the GTF file, to speed up processing
IndexGeneStartStop <- function(x){
    Counter <- 1
    gene.names <- NULL
    gene.starts <- NULL
    gene.stops <- NULL
    gene.starts[1] <- 1
    MAXrow <- length(x)
    for (i in seq(length(x))){
        # make a variable for the next gene, but make sure it isn't higher than the max 'i'
        nextValue <- min((i+1),MAXrow)
        if (x[i] != x[(i+1)] || i == length(x)){
            # annotate stop and fill in gene has table info
            gene.stops[Counter] <- i
            gene.names[Counter] <- x[i]
            # now determine if we continue (if we are not at the end of the GTF file)
            if (i != length(x)){
                Counter <- Counter + 1
                gene.starts[Counter] <- i+1
            }
        }
    }
    gene.data <- data.frame(gene.names,gene.starts,gene.stops)
    return(gene.data)
}

### This function transforms GFF annotation data tables to GRanges format
# We do not use rtracklayer because our GFF file contains 10 columns and may not be handled properly
GRangesFromGFF <- function(annot.GFF=NULL){
    # write out temporary GFF file
    write.table(annot.GFF,file=paste(tempdir(),"/temp.gtf3",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    # now import this file as a GFF3 file with rtracklayer
    GFF.GRanges <- rtracklayer::import(con=paste(tempdir(),"/temp.gtf3",sep=""),format = "GFF3")
    # clean up
    file.remove(paste(tempdir(),"/temp.gtf3",sep=""))
    # export this GRanges object & end function
    return(GFF.GRanges)
}

############################################################################
############################## process_Counts ##############################
############################################################################

### take rtracklayer objects & reformat them into standard GTF
GRangesToGFF <- function(GFF.GRanges=NULL){

    ### This function converts rtracklayer imported GTF files into flattened GFF files
    ### it functions by converting rtracklayer format back to GTF, and then calls GFF_convert()

    ### export rtracklayer object to GTF file (temporary gzip file)
    rtracklayer::export(GFF.GRanges, format = "GFF3", con=paste(tempdir(),"/temp.gff3",sep=""))

    ### now read this file back in
    annot.GFF <- read.table(file=paste(tempdir(),"/temp.gff3",sep=""),sep="\t", stringsAsFactors=FALSE)

    ### clean up
    file.remove(paste(tempdir(),"/temp.gff3",sep=""))

    ### return result & end function
    return(annot.GFF)
}

### reformat the GFF3 format to the format better suited for analysis
reformat_GFF3 <- function(annot.GFF=NULL){
    # make column 2 gene_id:exon_bin
    annot.GFF[,2] <- gsub(".*ID=(.*?);.*", "\\1", annot.GFF[,9])
    # make column 3 into name column
    annot.GFF[,3] <- gsub(".*Name=(.*?);.*", "\\1", annot.GFF[,9])
    # lastly make column 9 into transcripts column
    annot.GFF[,9] <- gsub(".*Transcripts=(.*?)", "\\1", annot.GFF[,9])
    # return & end function
    return(annot.GFF)
}

### normalize library size function (data must be normal space, not log2)
# any feature IDs (gene names, exon bins, etc.) must be in rownames, NOT data column
normalizeLibrarySizes <- function(raw.Data=NULL){
    # prevent factors
    options(stringsAsFactors=FALSE)
    # duplicate raw data original data
    original.Data <- raw.Data
    # add 1 count and take the log2 of reads
    raw.Data <- raw.Data+1
    raw.Data <- log2(raw.Data)
    # count columns and rows
    NCols <- ncol(raw.Data)
    NRows <- nrow(raw.Data)
    # remove genes/exons which have fewer than 8 reads (3 in log2 space) in any one condition
    raw.Data <- raw.Data[which(apply(raw.Data,1,min) >= 3),]

    # generate base mean
    raw.Data$baseMean <- apply(raw.Data,1,mean)
    # sizeFactors will be filled with % to adjust each column to baseMean
    mean.Diff=NULL
    # loop through conditions
    for (i in seq(NCols)){
        # difference between column expression & baseMean
        expr.Diff <- raw.Data[,i] - raw.Data$baseMean
        # now make a matrix to compute Lower & Upper percentages for shorth window
        Quant <- matrix(0,nrow=61,ncol=3)
        for (n in seq(61)){
            # Lower percentile of shorth window
            Lower <- 0.095 + n*0.005
            # Upper percentile of shorth window
            Upper <- Lower+0.5
            # Now compute quantile values for these percentile values & add them to Quant
            Quant[n,1:2] <- quantile(expr.Diff, c(Lower, Upper))
            # Compute the size of the shorth quantile window
            Quant[n,3] <- Quant[n,2] - Quant[n,1]
        }
        # Find the smallest shorth quantile window index
        minShorthIndex <- match(min(Quant[,3]),Quant[,3])
        # Now grab the shorth quantile data from that index
        Shorth <- Quant[minShorthIndex,c(1,2)]
        # Calculate the mean difference in log2 space between this raw.Data column & baseMean
        mean.Diff[i] <- 2^(mean(expr.Diff[which(expr.Diff >= Shorth[1] & expr.Diff <= Shorth[2])]))
    }
    ## now adjust raw.Data based on the sizeFactors we have generated
    adjusted.Data <- t(t(original.Data)/mean.Diff)
    # add row and column names
    rownames(adjusted.Data) <- rownames(original.Data)
    colnames(adjusted.Data) <- colnames(original.Data)[1:NCols]
    # now return output
    return(adjusted.Data)
}

parseAmbiguousReads <- function(read.Counts=NULL, annot.GFF=NULL){
    ### parse improper use of function
    if (is.null(annot.GFF) == TRUE){
        stop(call="You are attempting to run parseAmiguousReads outside of processCounts.
This function is designed to only be called internally by processCounts.
If you are running processCounts and still receiving this message, try re-generating your GFF file.")
    }

    ### generate GRanges
    adjcount.GRanges <- GRanges(seqnames=annot.GFF$V1, ranges=IRanges(as.numeric(annot.GFF$V4),
                                                                      as.numeric(annot.GFF$V5)))
    ### count overlaps per feature
    overlap.GRanges <- countOverlaps(adjcount.GRanges, adjcount.GRanges)

    ### set reads in overlapping exon bins from overlap.GRanges to 0, if any exist
    if (max(overlap.GRanges) > 1){
        read.Counts[which(overlap.GRanges > 1),] <- 0
    }
    return(read.Counts)
}


#####################################################################################################
#######################################  ExCluster Functions  #######################################
#####################################################################################################

#### Geometric mean function
gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#Compute effect size
Calc.effect.size <- function(Mean1,Mean2,Var1,Var2){
    Effect_size <- (Mean1-Mean2)/(sqrt((Var1+Var2)/2))
}

### function to determine number of clusters
Determine_HC_Distance <- function(x, HC_cutree){
    which(HC_cutree%in%x)
}

### using for loops, saving 1/2 time
Generate_ES_dists <- function(Mean1=NULL,Var1=NULL){
    x=length(Mean1)
    ES_mat <- matrix(0,nrow=x,ncol=x)

    # fill in nth row and jth column data for ES difference matrix
    for (n in seq(x-1)){
        for (j in seq((n+1),x)){
            ES_mat[n,j] <- Calc.effect.size(Mean1[n],Mean1[j],Var1[n],Var1[j])
            ES_mat[j,n] <- ES_mat[n,j] * -1
        }
    }
    return(ES_mat)
}

### function to compute the distances between all  clusters (apply max to its output)
Compare_HC_Distances <- function(NumClusters=NULL, ES_Dists=NULL, HC_indices=NULL){
    ES_Clust_Dist <- max(vapply(seq(NumClusters), FUN=function(z){
        obsDist <- mean(ES_Dists[HC_indices[[z]],-c(HC_indices[[z]])])
        return(obsDist)},FUN.VALUE = c(Res=0)))
    return(ES_Clust_Dist)
}


Permuted_ES_nullhypo <- function(r, Sim_log2FC=NULL, Sim_log2var=NULL, NumRows=NULL, DISTANCE=NULL, LINKAGE=NULL,
                                 NumClusters=NULL){

    pMean1 <- as.numeric(Sim_log2FC[r,])
    pVar1 <- as.numeric(Sim_log2var[r,])

    # ES difference matrix
    pES_permuted <- Generate_ES_dists(pMean1,pVar1)
    # ES distance matrix (Euclidean by default)
    pES_Dists <- as.matrix(dist(pES_permuted,method=DISTANCE))

    if (NumRows > 3){
        ##Generate hclust object
        HC <- hclust(as.dist(pES_Dists),method=LINKAGE)
        ##cut tree
        HC_cutree <- cutree(HC,k=NumClusters)
        ## find row indices per cluster
        HC_indices <- lapply(seq(NumClusters),FUN=function(y){which(HC_cutree==y)})
    }else{
        HC_indices <- seq(NumRows)
    }

    ### generate distance matrix from distance object
    pES_Clust_Dist <- max(vapply(seq(NumClusters), FUN=function(z){
        pDist <- mean(pES_Dists[HC_indices[[z]],-c(HC_indices[[z]])])
        return(pDist)},FUN.VALUE = c(Res=0)))
    return(pES_Clust_Dist)
}

### Note this function was adapted in full or in part from the NbClust R packaage
# NbClust offers no warranty on this function as per their GPL-2 license
centers <- function(cl, x) {
    x <- as.matrix(x)
    n <- length(cl)
    k <- max(cl)
    centers <- matrix(nrow = k, ncol = ncol(x))
    {
        for (i in seq(k)) {
            for (j in seq(ncol(x))) {
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
    for (u in seq(k)) cluster.size[u] <- sum(cl == u)
    for (u in seq(k)) {
        for (j in seq(ncol(x))) {
            for (i in seq(n)) {
                if (cl[i] == u)
                    variance.clusters[u, j] <- variance.clusters[u, j] + (x[i, j] - centers.matrix[u, j])^2
            }
        }
    }
    for (u in seq(k)) {
        for (j in seq(ncol(x))) variance.clusters[u, j] = variance.clusters[u, j]/cluster.size[u]
    }
    variance.matrix <- numeric(0)
    for (j in seq(ncol(x))) variance.matrix[j] = var(x[, j]) * (n - 1)/n
    Somme.variance.clusters <- 0
    for (u in seq(k)) Somme.variance.clusters <- Somme.variance.clusters +
        sqrt((variance.clusters[u, ] %*% (variance.clusters[u, ])))
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
    for (m in seq(MAX_cl)){
        N <- seq(MAX_cl)
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
Cutree_SD.Index <- function(DISTANCE=NULL, NRows=NULL, HC=NULL, ES_mat=NULL){
    ## hard coded at minimum of 2 clusters
    SD_values <- array()
    for(C in seq(2,NRows)){
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
    for (i in seq(28)){
        Lower <- (i)/400
        Upper <- (i+7)/400
        NumEstNull[i+1] <- length(x[which(x >= Lower & x < Upper)])
    }
    # check to see if p-values have hit a low plateau (null hypothesis)
    for (i in seq(20)){
        # Check maximum number of pvals in the 4 below the current i
        NumBelow <- NumEstNull[seq((i+1),(i+5))]
        # Convert these to TRUE/FALSE
        NumBelowBoolean <- NumEstNull[i] > NumBelow
        # Now see how many of these bins the current NumEstNull[i] is less than
        LessThan <- length(which(NumBelowBoolean == FALSE))
        if (LessThan >= 2) break
    }
    NumEstNull <- min(i/400,0.05)
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
    for (n in seq(length(x))){
        # check to make sure that we aren't on the first value, and the previous value isn't FDR == 1
        if (n > 1 && FDRvalues[(n-1)] >= 0.99) {
            FDRvalues[n] <- 1
        }else{
            ### now generate a null hypothesis 2 times to estimate the NumExp and Num Obs
            NumExp <- NULL
            NumObs <- NULL
            for (h in seq(3)){
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
    for (n in seq(2,FDRlim[1])){
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
        MinFDR <- min(FDRvalues[seq((n-1),MaxN)])
        # now replace the FDRvalues[n] if FDRvalues[n] > Min FDR
        if (FDRvalues[(n-1)] > MinFDR){
            FDRvalues[(n-1)] <- MinFDR
        }
    }
    return(FDRvalues)
}

