processCounts <- function(bam.Files=NULL,sample.Names=NULL,annot.GFF=NULL,GFF.File=NULL,pairedReads=FALSE,out.File=NULL){

    ### make sure that bam.Files were entered when running processCounts
    if(is.null(bam.Files) == TRUE){
        stop(call="No BAM files were entered when running processCounts(). Please ensure the bam.Files argument is assigned a character vector of full BAM file paths when calling processCounts().
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

    # check to make sure either the GFF data is present, or a GFF file path is given.
    if(is.null(annot.GFF) == TRUE){
        if (is.null(GFF.File) == FALSE){
            if (file.exists(GFF.File) == TRUE){
                annot.GFF <- read.table(file=GFF_file,header=FALSE,stringsAsFactors=FALSE)
            }else{
                stop(call="The GFF file path you provided to GFF.File does not exist. Please verify you have given the correct, full filepath including file extension. For example, a file path may look like: /Users/username/path/to/file.gff")
            }
        }else{
            stop(call="The GFF.File argument was not provided when processCounts was run, nor was the GFF.Data argument assigned. One of these two arguments must be assigned.
                The GFF.File argument should be a full file path to a GFF file, such as: /Users/username/path/to/file.gff
                Alternatively, you can give the annot.GFF argument a GFF data frame from the GFF_convert function, such as: annot.GFF = GFF")
        }
    }

    ### now load Rsubread library if it isn't loaded
    if (("Rsubread" %in% rownames(installed.packages())) == FALSE){
        library(Rsubread)
        # now make sure it loaded again -- if not, error out this function
        if (("Rsubread" %in% rownames(installed.packages())) == FALSE){
            stop(call="The Rsubread package from Bioconductor is required to run this function.
It appears that attempts to load this package have failed. Please re-install Rsubread from Bioconductor.
If you are unable to manually install & load Bioconductor or Rsubread, please contact the appropriate maintainer.")
        }
    }


    # geometric mean function
    gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }

    ### Convert GFF annotations to SAF
    SAF.annot <- data.frame(annot.GFF$V2,annot.GFF$V1,annot.GFF$V4,annot.GFF$V5,annot.GFF$V7)
    colnames(SAF.annot) <- c("GeneID","Chr","Start","End","Strand")

    ### Run featureCounts on BAM files
    ### The authors and original license holders of featureCounts and the Rsubread package make no warranty for its performance
    fC <- featureCounts(files = bam.Files, annot.ext = SAF.annot, isGTFAnnotationFile = FALSE, requireBothEndsMapped = FALSE,
                        allowMultiOverlap = TRUE, strandSpecific = 1, largestOverlap = FALSE, isPairedEnd = pairedReads)

    cat('',"Running library size normalization...",'',sep="\n")

    ### clean up SAF.annot (not needed anymore)
    rm(SAF.annot)

    DataTable <- data.frame(fC$counts)
    colnames(DataTable) <- c(sample.Names)

    NRows <- nrow(DataTable)
    NCols <- ncol(DataTable)

    ### Save original table
    OriginalTable <- DataTable

    ### add 1 count to each exon and take the log2, also calculate the geometric baseMean
    DataTable <- DataTable+1
    DataTable <- log2(DataTable[,1:NCols])
    DataTable$baseMean <- apply(DataTable[,1:NCols],1,gm_mean)

    #### Remove exons which are not well expressed
    DataTable <- DataTable[which(DataTable$baseMean >= 4),]
    NRows <- nrow(DataTable)

    Indices <- matrix(data=0,nrow=NRows,ncol=NCols)
    rownames(Indices) <- rownames(DataTable)

    ### Make sure we have at least 2000 expressed exon bins to normalize
    if (NRows > 2000){
        for (i in 1:NCols){
            #working <- get(Condition[i])
            a <- rep(0,NRows)
            working <- cbind(DataTable,a)
            working <- cbind(DataTable,foldChange=a)
            working$foldChange <- DataTable[,i] - DataTable$baseMean
            Quant <- matrix(0,nrow=41,ncol=3)
            for (n in 1:41){
                Lower <- 0.09 + n*0.01
                Upper <- Lower+0.4
                Quant[n,1:2] <- quantile(working$foldChange, c(Lower, Upper))
                Quant[n,3] <- Quant[n,2] - Quant[n,1]
            }
            Minimum <- match(min(Quant[,3]),Quant[,3])
            Shorth <- Quant[Minimum,]
            Indices[which(working$foldChange >= Shorth[1] & working$foldChange <= Shorth[2]),i] <- 1

            rm(working)
            rm(a)
            rm(Quant)
        }

        ### now add up the values in the 'Indices' dataframe, where 1 = within the shorth for that sample
        ### rows with a sum of 6 had that exon within the shorth for every exon, indicating they are good 'control' exons
        indicesSum <- apply(Indices,1,sum)
        CommonData <- DataTable[which(indicesSum == 6),]

        sizeFactors <- matrix(nrow=1,ncol=NCols)

        ### now restore original read counts to be adjusted in the new 'AdjustedData'
        AdjustedData <- OriginalTable
        rm(OriginalTable)
        rm(DataTable)

        for (i in 1:NCols){
            templog2FC <- CommonData[,i] - CommonData$baseMean
            sizeFactors[1,i] <- 2^(median(templog2FC))
            AdjustedData[,i] <- AdjustedData[,i]/sizeFactors[1,i]
        }

        ### if we don't have > 1000 expressed exon bins, run simple normalization
    }else{
        # clean up
        rm(DataTable)

        # simple library size normalization
        sampleSums <- apply(OriginalTable,2,sum)
        sampleBaseMean <- mean(sampleSums)
        sizeFactors <- sampleSums/sampleBaseMean

        # make results AdjustedData table
        AdjustedData <- OriginalTable
        # clean up
        rm(OriginalTable)
        # divide reads by library sizeFactors
        AdjustedData <- t(t(AdjustedData)/sizeFactors)

        # warn the user that they didn't have a good amount of data
        warning(call="You had fewer than 1000 exon bins with RNA-seq reads counted into them.
Please be aware that ExCluster expects whole transcriptome GFF and BAM files counted, or else it cannot properly correct library sizes and tune statistics.
If you are simply testing an ExCluster example, please ignore this warning.")
    }

    if (is.null(out.File) == FALSE){
         write.table(AdjustedData[,1:NCols],file=out.File,row.names=FALSE,col.names=TRUE,quote=FALSE)
    }
    return(AdjustedData)

    cat('',"processCounts function has completed.",sep="\n")
}
