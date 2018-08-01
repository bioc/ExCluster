GRangesFromGFF <- function(annot.GFF=NULL, GFF.File=NULL){

    ############################################################################################
    ############################## Very GFF file format integrity ##############################
    ############################################################################################

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
Please ensure you are using a GFF file from GFF_convert, and any saved GFF file is read in without row/column names.
See the ExCluster manual and vignette for more information.")
        }
    }

    ### check to make sure that each element of column 3 in the GFF file begins with 'exonic_part'
    exonic_part.Values <- unique(substr(annot.GFF[,3],1,11))
    # if the length of this unique value is > 1, or the result does not == 'exonic_part', throw an error
    if (length(exonic_part.Values) > 1){
        stop(call = "Your GFF file did not have column 3 beginning with 'exonic_part' in all cases.
Please ensure your GFF file is produced from the GFF_convert function, and read into R without row/column names.
If this error persists, please re-produce your GFF file from a GTF file using GFF_convert.")
    }

    ### check to make sure columns 4 & 5 are numeric start/stops with positive differences only
    # we do this by subtracting 'start' from 'stop' and making sure the math works, and minimum = 0
    # test to see if the math will work (will only work if columns 4 & 5 are properly formatted)
    coord.Diff <- try(as.numeric(annot.GFF[,5]) - as.numeric(annot.GFF[,4]),silent = TRUE)
    # check to see if we have an error
    if (substr(coord.Diff[1],1,5) == "Error"){
        stop(call = "The start and stop columns (4 & 5) of your GFF file were not properly formatted.
Please ensure your GFF annotations were produced from the GFF_convert function.
If you are reading saved GFF annotations from a file, please ensure no row or column names were applied to the file.
If this error persists, please re-produce your GFF file from a GTF file using GFF_convert.")
        }
    # if that did not fail, the code continues and we verify that all stop - start subtractions are >= 0
    if (min(coord.Diff) < 0){
        stop(call = "There were negative distances between your stop and start locations in your GFF file.
Please ensure your GFF annotations were produced from the GFF_convert function.
If you are reading saved GFF annotations from a file, please ensure no row or column names were applied to the file.
Columns 4 & 5 of the GFF file format should contain start/stop locations with only numeric values throughout.
If this error persists, please re-produce your GFF file from a GTF file using GFF_convert.")
    }

    ### check to make sure column 7 unique values have at least "+" or "-" (could be just one for test GFF files)
    # unique values of column 7 in the GFF file
    strand.Values <- unique(annot.GFF[,7])
    # is "+" present? if true this will == 1, if false this will == 0
    positive.Length <- length(which('+'%in%strand.Values))
    # is "-" present? if true this will == 1, if false this will == 0
    negative.Length <- length(which('-'%in%strand.Values))
    # make sure there are no more than 2 unique values of column 7
    if (length(strand.Values) > 2){
        stop(call = "The strand column of your GFF file (column 7) had more than 2 unique values.
This column should only have '+' and '-' values. Please ensure your GFF file is read into R without row/column names.
If this error persists, please re-produce your GFF file from a GTF file using GFF_convert.")
    }
    # if the previous check passed, make sure at least one of '-' or '+' is present
    if ((positive.Length + negative.Length) < 1){
        stop(call = "The strand column of your GTF file (column 7) did not contain '+' or '-' values.
This column should only have '+' and '-' values. Please ensure your GFF file is read into R without row/column names.
If this error persists, please re-produce your GFF file from a GTF file using GFF_convert.")
    }

    #############################################################################################
    ############################# Now convert GFF format to GRanges #############################
    #############################################################################################

    # convert to GRanges
    GFF.GRanges <- GRanges(seqnames=annot.GFF$V1, ranges=IRanges(as.numeric(annot.GFF$V4),as.numeric(annot.GFF$V5)),
                           exon_bin_id=annot.GFF$V2, exonic_part=annot.GFF$V3, transcript_ids=annot.GFF$V9,
                           gene_name=annot.GFF$V10)

    # export this GRanges object & end function
    return(GFF.GRanges)
}
