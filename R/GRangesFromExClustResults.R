GRangesFromExClustResults <- function(results.ExClust=NULL){

    ############################################################################################
    ############################## Very ExClust Results integrity ##############################
    ############################################################################################

    ### make sure column names are correct (all columns up until read count data)
    expected.Names <- c("EnsG","Exon_bin","Gene_name","Cluster","log2FC","log2Var","pval","FDR",
                        "Chr","Strand","Start","End")

    ### make sure that results.ExClust is a data frame with at least 2+ columns
    if (length(ncol(results.ExClust)) != 1){
        if(ncol(results.ExClust) < 2){
            stop(call="You did not enter a valid data structure for results.ExClust with columns.
Please ensure you pass results.ExClust a data structure from the main ExCluster function.
This data structure should have 12 columns at minimum, plus more if read count columns are present.")
        }
    }

    ### make sure that results.ExClust is a data frame with at least 2+ rows
    if (length(nrow(results.ExClust)) != 1){
        if(nrow(results.ExClust) < 2){
            stop(call="You did not enter a valid data structure for results.ExClust with more than 1 row.
Please ensure you pass results.ExClust a data structure from the main ExCluster function.
This data structure should have thousands of rows (exon bins) for any organism transcriptome.")
        }
    }

    ### now check to see if these expected names match actual names
    equal.Vectors <- all.equal(colnames(results.ExClust)[seq(length(expected.Names))], expected.Names)
    # if these vectors are not equal, throw an error
    if (equal.Vectors == FALSE){
        stop(call="The first 12 column names of results.ExClust did not match the expected format.
Please ensure you are passing the ExCluster main function results table directly to results.ExClust.")
    }

    ### check to make sure the Start/End columns are numeric with positive differences only
    # test to see if the math will work (will not work if columns are not numeric)
    coord.Diff <- try(as.numeric(results.ExClust$End) - as.numeric(results.ExClust$Start),silent = TRUE)
    # check to see if we have an error
    if (substr(coord.Diff[1],1,5) == "Error"){
        stop(call = "The Start/End columns of your ExClust results table were not properly formatted.
These columns should consist of only numeric values, but instead contained non-numeric values.
Please ensure you are passing the ExCluster main function results table directly to results.ExClust.")
        }
    # if that did not fail, the code continues and we verify that all stop - start subtractions are >= 0
    if (min(coord.Diff) < 0){
        stop(call = "There were negative distances between your Start/End locations in your ExClust results.
These distances should always be positive when subtracting Start locations from End locations, per exon bin.
Please ensure you are passing the ExCluster main function results table directly to results.ExClust.")
    }

    ### check to make sure the Strand column unique values have at least "+" or "-"
    # unique values of Strand column in the ExClust results
    strand.Values <- unique(results.ExClust$Strand)
    # is "+" present? if true this will == 1, if false this will == 0
    positive.Length <- length(which('+'%in%strand.Values))
    # is "-" present? if true this will == 1, if false this will == 0
    negative.Length <- length(which('-'%in%strand.Values))
    # make sure there are no more than 2 unique values contained within the Strand column
    if (length(strand.Values) > 2){
        stop(call = "The Strand column of your ExClust results contained more than 2 unique values.
This column should only have '+' and '-' values.
Please ensure you are passing the ExCluster main function results table directly to results.ExClust.")
    }
    # if the previous check passed, make sure at least one of '-' or '+' is present
    if ((positive.Length + negative.Length) < 1){
        stop(call = "The Strand column of your ExClust results did not contain '+' or '-' values.
This column should only have '+' and '-' values.
Please ensure you are passing the ExCluster main function results table directly to results.ExClust.")
    }

    #############################################################################################
    ########################## Now convert ExClust Results to GRanges ###########################
    #############################################################################################

    # convert to GRanges
    ExClustResults.GRanges <- GRanges(seqnames=results.ExClust$Chr,
                                      ranges=IRanges(as.numeric(results.ExClust$Start),as.numeric(results.ExClust$End)),
                                      EnsG=results.ExClust$EnsG, Exon_bin=results.ExClust$Exon_bin,
                                      Gene_name=results.ExClust$Gene_name, Cluster=results.ExClust$Cluster,
                                      log2FC=results.ExClust$log2FC, log2Var=results.ExClust$log2Var,
                                      pval=results.ExClust$pval, FDR=results.ExClust$FDR,
                                      counts=results.ExClust[,c(13:ncol(results.ExClust))])

    # export this GRanges object & end function
    return(ExClustResults.GRanges)
}
