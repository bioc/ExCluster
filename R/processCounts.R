processCounts <- function(bam.Files=NULL,sample.Names=NULL,annot.GFF=NULL,GFF.File=NULL,pairedReads=FALSE,out.File=NULL,  temp.Dir=NULL, num.Cores=NULL){

    ### make sure R doesn't add factors -- R creators never should have made this default = TRUE
    options(stringsAsFactors=FALSE)

    ### make sure that bam.Files were entered when running processCounts
    if(is.null(bam.Files) == TRUE){
        stop(call="No BAM files were entered when running processCounts().
Please ensure the bam.Files argument is assigned a character vector of full BAM file paths when calling processCounts().
Alternatively, you can assign the BAM file names as a character vector to the bam.Files variable, such as:
    bam.Files <- c('/path/to/file1.bam','/path/to/file2.bam',...)")
    }

    ### now make sure the sample.Names variable is set
    if(is.null(sample.Names) == TRUE){
        stop(call="No sample names were entered to correspond to your BAM files.
             Please make sure the sample.Names argument is assigned a character vector when calling processCounts().
             Alternatively, you can assign the character vector to the sample.Names variable as follows:
             sample.Names <- c('control1_rep1','control1_rep2',...)")
    }

    ### make sure that 'pairedReads' variable is logical
    if(is.logical(pairedReads) == FALSE){
        stop(call="You did not correctly enter the 'pairedReads' variable as TRUE or FALSE.
             Please use pairedReads=TRUE or pairedReads=FALSE when calling processCounts()")
    }

    ### Check to make sure the number of sample & sample names are the same
    if (length(bam.Files) != length(sample.Names)){
        stop(call="The number of BAM files and sample names are not the same!
             Please ensure the lengths of bam.Files and sample.Names variables are the same.")
    }


    ########################################## Library size Normalization Function ##########################################

    ## columns must be "conditions" and rows must be features, i.e. genes or exons
    # Read counts must be in normal space (i.e. not log2)
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
        # remove genes/exons which have fewer than 8 reads in any one condition
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
            Quant <- matrix(0,nrow=31,ncol=3)
            for (n in seq(31)){
                # Lower percentile of shorth window
                Lower <- 0.09 + n*0.01
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
            Shorth <- Quant[minShorthIndex,c(1:2)]
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

    ########################################## Read counts with Rsubread ##########################################

    # check to make sure either the GFF data is present, or a GFF file path is given.
    if(is.null(annot.GFF) == TRUE){
        if (is.null(GFF.File) == FALSE){
            if (file.exists(GFF.File) == TRUE){
                annot.GFF <- read.table(file=GFF.File,header=FALSE,stringsAsFactors=FALSE)
            }else{
                stop(call="The GFF file path you provided to GFF.File does not exist.
Please verify you have given the correct, full filepath including file extension.
    For example, a file path may look like: /Users/username/path/to/file.gff")
            }
        }else{
            stop(call="The GFF.File argument was not provided when processCounts was run, nor was the annot.GFF argument assigned.
One of annot.GFF or GFF.File must be specified when running processCounts().
    The GFF.File argument should be a full file path to a GFF file, such as: /Users/username/path/to/file.gff
    Alternatively, you can give the annot.GFF argument a GFF data frame from the GFF_convert function, such as: annot.GFF = GFF")
        }
    }
    ### Convert GFF annotations to SAF
    SAF.annot <- data.frame(annot.GFF$V2,annot.GFF$V1,annot.GFF$V4,annot.GFF$V5,annot.GFF$V7)
    colnames(SAF.annot) <- c("GeneID","Chr","Start","End","Strand")

    ### default nCores = 1 if not specified, or incorrectly specified
    nCores <- 1
    # now check if we are using multiple cores
    if (is.null(num.Cores) == FALSE){
        # check to make sure num.Cores is a 'double' type variable
        if (is.double(num.Cores) == TRUE){
            nCores <- num.Cores
        }else{
            warning(call="You entered a non-integer value for num.Cores -- it should have a value of >= 1.
    Defaulting number of cores to 1 for this run.")
            nCores <- 1
        }
    }

    ### default the 'tmpDir' argument for featureCounts to tempdir() if not specified
    tmpDir <- tempdir()
    # now check if the temp.Dir argument was specified in ExClusters input
    if (is.null(temp.Dir) == FALSE){
        # make sure we can write to temp.Dir
        WriteCheck <- file.access(temp.Dir, mode=2)
        if (WriteCheck != 0){
            tmpDir <- temp.Dir
        }else{
            warning(call = "The temp.Dir path you provided was not accessible by R for read/write.
    Please ensure the temp.Dir argument is given a valid folder path with read/write permissions.
    Defaulting temporary directory to R's tempdir() function.")
        }
    }

    ### Run featureCounts on BAM files
    ### The authors and original license holders of featureCounts and the Rsubread package make no warranty for its performance
    fC <- featureCounts(files = bam.Files, annot.ext = SAF.annot, isGTFAnnotationFile = FALSE, nthreads = nCores, tmpDir=tmpDir,
            requireBothEndsMapped = FALSE, allowMultiOverlap = TRUE, largestOverlap = FALSE, isPairedEnd = pairedReads)

    cat('',"Running library size normalization...",'',sep="\n")

    ### clean up SAF.annot (not needed anymore)
    rm(SAF.annot)

    DataTable <- data.frame(fC$counts)
    colnames(DataTable) <- c(sample.Names)

    ########################################## Now normalize read counts ##########################################

    # check if we have at least 1000 rows of data to normalize counts to
    if (nrow(DataTable[which(apply(DataTable,1,min) >= 8),]) < 1000){
        # warn the user that they didn't have a good amount of data if less than 1000 elibile features (i.e. exon bins)
        warning(call="You had fewer than 1000 exon bins with RNA-seq reads counted into them across conditions.
Please be aware that ExCluster expects whole transcriptome GFF and BAM files counted, not individual chromosome data.
Without these full files, ExCluster cannot properly correct library sizes and tune statistics.
    >> If you are simply testing an ExCluster example, please ignore this warning.")

        # simple library size normalization
        sampleSums <- apply(DataTable,2,sum)
        # compute baseMean library size (across all conditions)
        sampleBaseMean <- mean(sampleSums)
        # compute sizeFactor differences for each sample vs basemean
        sizeFactors <- sampleSums/sampleBaseMean
        # now adjust DataTable based on these size factors
        adjusted.Counts <- t(t(DataTable)/sizeFactors)

    }else{
        adjusted.Counts <- normalizeLibrarySizes(DataTable)
    }

    ########################################## Final processing & output  ##########################################

    if (is.null(out.File) == FALSE){
         write.table(adjusted.Counts,file=out.File,sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
    }
    return(adjusted.Counts)

    cat('',"processCounts function has completed.",sep="\n")
}
