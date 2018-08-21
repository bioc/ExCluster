processCounts <- function(bam.Files=NULL,sample.Names=NULL,annot.GFF=NULL,GFF.File=NULL, pairedReads=FALSE,
                          stranded.Reads=FALSE, out.File=NULL, temp.Dir=NULL, num.Cores=NULL){

    ### make sure R doesn't add factors -- R creators never should have made this default = TRUE
    options(stringsAsFactors=FALSE)
    # make sure that bam.Files were entered when running processCounts
    if(is.null(bam.Files) == TRUE){
        stop(call = ExCluster_errors$BAM_argument_missing)
    }
    #now make sure the sample.Names variable is set
    if(is.null(sample.Names) == TRUE){
        stop(call = ExCluster_errors$sample_name_argument_missing)
    }
    # make sure that 'pairedReads' variable is logical
    if(is.logical(pairedReads) == FALSE){
        stop(call = ExCluster_errors$paired_reads_not_logical)
    }
    # Check to make sure the number of sample & sample names are the same
    if (length(bam.Files) != length(sample.Names)){
        stop(call = ExCluster_errors$BAM_sample_length_unequal)
    }

    ########################################## Read counts with Rsubread ##########################################

    ### check to make sure either the GFF data is present, or a GFF file path is given.
    if(is.null(annot.GFF) == TRUE){
        if (is.null(GFF.File) == FALSE){
            if (file.exists(GFF.File) == TRUE){
                annot.GFF <- read.table(file=GFF.File,header=FALSE,stringsAsFactors=FALSE)
            }else{
                stop(call = ExCluster_errors$GFF_file_doesnt_exist)
            }
        }else{
            stop(call = ExCluster_errors$GFF_file_argument_missing)
        }
    }else{
        ### now check to see if annot.GFF is an S4 object (GRanges) & reformat it
        if (substr(typeof(annot.GFF),1,2) == "S4"){
            annot.GFF <- GRangesToGFF(annot.GFF)
        }
    }

    ### check to make sure the GFF file has 9 columns
    if (ncol(annot.GFF) != 9){
        stop(call=ExCluster_errors$GFF_missing_columns)
    }

    ### now check to make sure that column 9 contains 'ID=', 'Name=', and 'Transcripts=' strings
    ID.indices <- length(grep(pattern = "ID=", as.character(annot.GFF[seq(2),9])))
    Name.indices <- length(grep(pattern = "Name=", as.character(annot.GFF[seq(2),9])))
    Transcripts.indices <- length(grep(pattern = "Transcripts=", as.character(annot.GFF[seq(2),9])))
    # if the lengths of these indices do not sum to 6, throw an error
    if (ID.indices+Name.indices+Transcripts.indices != 6){
        stop(call=ExCluster_errors$GFF_missing_GFF3_fields)
    }

    ### now reformat the GFF file
    annot.GFF <- reformat_GFF3(annot.GFF)

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
            warning(call=ExCluster_errors$non_integer_numcores)
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
            warning(call = ExCluster_errors$temp_dir_inaccessible)
        }
    }

    ### are reads stranded?
    if (stranded.Reads == FALSE){
        stranded.Reads <- 0
    }else{
        stranded.Reads <- 1
    }

    ### Run featureCounts on BAM files
    ### The authors and original license holders of featureCounts and the Rsubread package make no warranty for its performance
    fC <- featureCounts(files = bam.Files, annot.ext = SAF.annot, isGTFAnnotationFile = FALSE, nthreads = nCores, tmpDir=tmpDir,
            requireBothEndsMapped = FALSE, allowMultiOverlap = TRUE, largestOverlap = FALSE, isPairedEnd = pairedReads,
            strandSpecific = stranded.Reads)

    ### run ambiguous read removal if stranded.Reads == FALSE
    if (stranded.Reads == 0){
        # set ambiguous reads to zero if reads are unstranded
        fC <- parseAmbiguousReads(read.Counts = fC, annot.GFF=annot.GFF)
    }

    cat('',"Running library size normalization...",'',sep="\n")

    ### clean up SAF.annot (not needed anymore)
    rm(SAF.annot)

    DataTable <- data.frame(fC$counts)
    colnames(DataTable) <- c(sample.Names)

    ########################################## Now normalize read counts ##########################################

    # check if we have at least 1000 rows of data to normalize counts to
    if (nrow(DataTable[which(apply(DataTable,1,min) >= 8),]) < 1000){
        # warn the user that they didn't have a good amount of data if less than 1000 elibile features (i.e. exon bins)
        warning(call = ExCluster_errors$low_exon_bin_count)

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

    ### write output file if necessary
    if (is.null(out.File) == FALSE){
         write.table(adjusted.Counts,file=out.File,sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
    }

    ### return counts and finish function
    return(adjusted.Counts)
    cat('',"processCounts function has completed.",sep="\n")
}
