ExCluster <- function(exonCounts=NULL, cond.Nums=NULL, annot.GFF = NULL, GFF.File=NULL, out.Dir=NULL,
                      result.Filename=NULL, CombineExons=TRUE,plot.Results=FALSE,FDR.cutoff=0.05){

    ### make sure R doesn't add factors -- R creators never should have made this default = TRUE
    options(stringsAsFactors=FALSE)

    ### ExCluster progression message
    message("","Initializing ...",sep="\n")

    ### check to make sure a normalized exon count matrix was provided as input
    if (is.null(exonCounts) == TRUE){
        stop(call = ExCluster_errors$exon_counts_missing)
    }

    ### check to make sure condition numbers were specified
    if (is.null(cond.Nums) == TRUE){
        stop(call = ExCluster_errors$cond_nums_missing)
    }

    ### ensure the number of columns in the exonCounts parameter equal the number of cond.Nums
    if (ncol(exonCounts) != length(cond.Nums)){
        stop(call = ExCluster_errors$exon_count_length_incorrect)
    }

    ## grabbing unique condition numbers (should only ever be two conditions)
    uCondNums <- unique(cond.Nums)

    ## ensure exactly 2 conditions
    if (length(uCondNums) != 2){
        stop(call = ExCluster_errors$improper_cond_nums)
    }

    ## condition indices
    Indices1 <- grep(uCondNums[1],cond.Nums)
    Indices2 <- grep(uCondNums[2],cond.Nums)

    ## Ensure that both conditions have at least 2 replicates provided, or ExCluster will not function.
    if (length(Indices1) == 1){
        stop(call = ExCluster_errors$not_enough_replicates_cond1)
    }
    if (length(Indices2) == 1){
        stop(call = ExCluster_errors$not_enough_replicates_cond2)
    }

    ### check to ensure that there are no duplicated data columns in exonCounts
    CheckDups <- duplicated(t(exonCounts))
    if (any(CheckDups) == TRUE){
        DupCols <- which(CheckDups == TRUE)
        message(paste("Columns",DupCols,"were duplicated in your normalized exon count data.
    Please re-check your BAM file names and re-run processCounts."))
        stop(call = ExCluster_errors$exon_read_counts_duplicated)
        }

    ### check to make sure either a GFF data structure or GFF file path were specified
    if (is.null(annot.GFF) == TRUE){
        if (is.null(GFF.File) == TRUE){
            stop(call = ExCluster_errors$GFF_annotations_missing)
        }else{
            ### Test if GFF.File works
            if(file.exists(GFF.File) != TRUE){
                stop(call = ExCluster_errors$GFF_file_inaccessible)
            }
        }
    }

    #####################################################################################################
    #################################  READ COUNT & DATA PREPARATION  ###################################
    #####################################################################################################

    ### ExCluster progression message
    message('','Reading in data ...','',sep="\n")

    log2FC <- data.frame(exonCounts)
    rm(exonCounts)

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
        ### now check to see if annot.GFF is an S4 object (GRanges) & reformat it
        if (substr(typeof(annot.GFF),1,2) == "S4"){
            annot.GFF <- GRangesToGFF(annot.GFF)
        }
    }

    ######################### NOW VERIFY GFF FILE FORMAT INTEGRITY #########################

    if (ncol(annot.GFF) != 9){
        ### inform the user that the GFF file was not properly formatted
        stop(call= ExCluster_errors$ExClust_GFF_bad_colnum)
    }

    ### check to make sure that each element of column 3 in the GFF file begins with 'exon'
    exon_feature.Values <- unique(annot.GFF[,3])
    # if the length of this unique value is > 1, or the result does not == 'exon', throw an error
    if (length(exon_feature.Values) > 1){
        stop(call = ExCluster_errors$GFF_missing_exon_field)
    }

    ### check to make sure columns 4 & 5 are numeric start/stops with positive differences only
    # we do this by subtracting 'start' from 'stop' and making sure the math works, and minimum = 0
    # test to see if the math will work (will only work if columns 4 & 5 are properly formatted)
    coord.Diff <- try(as.numeric(annot.GFF[,5]) - as.numeric(annot.GFF[,4]),silent = TRUE)
    # check to see if we have an error
    if (substr(coord.Diff[1],1,5) == "Error"){
        stop(call = ExCluster_errors$GFF_numeric_error)
        }
    # if that did not fail, the code continues and we verify that all stop - start subtractions are >= 0
    if (min(coord.Diff) < 0){
        stop(call = ExCluster_errors$GFF_negative_dists)
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
        stop(call = ExCluster_errors$bad_GFF_strand_column)
    }
    # if the previous check passed, make sure at least one of '-' or '+' is present
    if ((positive.Length + negative.Length) < 1){
        stop(call = ExCluster_errors$GFF_plus_minus_missing)
    }
    # now check to make sure that column 9 contains 'ID=', 'Name=', and 'Transcripts=' strings
    ID.indices <- length(grep(pattern = "ID=", as.character(annot.GFF[seq(2),9])))
    Name.indices <- length(grep(pattern = "Name=", as.character(annot.GFF[seq(2),9])))
    Transcripts.indices <- length(grep(pattern = "Transcripts=", as.character(annot.GFF[seq(2),9])))
    # if the lengths of these indices do not sum to 6, throw an error
    if (ID.indices+Name.indices+Transcripts.indices != 6){
        stop(call=ExCluster_errors$GFF_missing_GFF3_fields)
    }

    ### if all of these checks passed, reformat the GFF3 formatted file into the prefered GFF format
    annot.GFF <- reformat_GFF3(annot.GFF)

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
    log2FC$Genes <- annot.GFF$V3
    log2FC$Bins <- log2Bins

    # compute means and variances for the log2FC data frame
    log2FC$log2FC <- rowMeans2(as.matrix(log2(log2FC[,Indices2]+1))) -
        rowMeans2(as.matrix(log2(log2FC[,Indices1]+1)))
    log2FC$log2Variance <- rowVars(as.matrix(log2(log2FC[,Indices2]+1))) +
        rowVars(as.matrix(log2(log2FC[,Indices1]+1)))

    ### now we work on organizing the newlog2FC read count data
    # remove all exon bin rows with fewer than 2 max reads (log2 reads >= 1) across samples
    newlog2FC <- newlog2FC[which(rowMaxs(as.matrix(newlog2FC)) >= 1),]
    # now add EnsID to remove duplicated rows
    newlog2FC$EnsID <- sub(":.*","",rownames(newlog2FC))
    # identify duplicated rows by condition means
    dupRows <- which(duplicated(newlog2FC))

    # remove duplicated rows if there are any
    if (length(dupRows) > 0){
        newlog2FC <- newlog2FC[-c(dupRows),]
    }
    # remove EnsiD column
    newlog2FC <- newlog2FC[,-c(ncol(newlog2FC))]
    # set up newlog2FC "Bins"
    newBins <- sub(".*:","",rownames(newlog2FC))

    # take Ensembl ID numbers for easy parsing of newlog2FC data
    newIDs <- parseGeneIDs(gene.Annot =  rownames(newlog2FC))

    ### Add Ensembl ID column to annot.GFF data frame for later sorting
    annot.GFF$EnsID <- parseGeneIDs(gene.Annot = annot.GFF[,2])

    ### remove ENSG00000 from ENSEMBL IDs to make them easier to sort (for looping through each unique ID)
    IDs <- parseGeneIDs(gene.Annot = rownames(log2FC))

    ### Check to make that annot.GFF EnsID column and IDs vector are 100% identical
    EqualVectors <- all.equal(annot.GFF$EnsID,IDs)
    if (EqualVectors != TRUE){
        stop(call=ExCluster_errors$match_GFF_counts_failed)
    }

    ###############################  IF TRUE, COMBINE COMMON EXONS  #####################################

    if (CombineExons == TRUE){
        ### ExCluster progression message
        message("Combining common exons ...","",sep="\n")

        # Count number of rows in annot.GFF
        NROW_GFF <- nrow(annot.GFF)
        # Counter for counting new summed exon bins
        Counter <- 1
        # Start point for first newIDs gene row number
        newIDstart <- 1
        # Variable start will be used through loop
        Start <- 1
        # make newlog2FC list to dump results from CombineExons into
        temp_log2FC <- vector("list", 1)

        for (i in seq(NROW_GFF)){
            nextValue <- min((i+1), NROW_GFF)
            ### check if we are on a new gene
            if ((IDs[i] != IDs[nextValue]) == TRUE | i == NROW_GFF){

                # stop is the current ith trow
                Stop <- i
                # Rows of the log2FC data frame
                log2Rows <- seq(Start,Stop)

                ### now determine if the new gene is in newIDs
                if ((IDs[Start] == newIDs[newIDstart]) == TRUE){
                    # temporary newlog2FC data
                    TempIDs <- newIDs[seq(newIDstart,(newIDstart+length(log2Rows)))]
                    # determine the Stop point of newIDs
                    newIDstop <- newIDstart+length(which(TempIDs%in%IDs[Start]))-1
                    # Rows of the newIDs vector
                    newlog2Rows <- seq(newIDstart,newIDstop)
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
                        for (j in seq(length(uExTrans))){
                            uIndices <- which(ExTrans%in%uExTrans[j])
                            if (length(uIndices) > 1){
                                TransMatrix[j,] <- colSums2(as.matrix(TempReads[uIndices,]))
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
                if (i != NROW_GFF){
                    # now start the next gene
                    Start <- i+1
                }
            }
        }
        newlog2FC <- do.call(rbind,temp_log2FC)
        # clean up
        rm(temp_log2FC)
    }
    # clean up annot.GFF
    rm(annot.GFF)
    rm(log2Bins)
    rm(IDs)
    rm(newIDs)

    ################################### Now run ExCluster main function ##################################

    log2Clusters <- ExClust_main_function(newlog2FC, Indices1, Indices2, NumReps1, NumReps2)

    #################################### Compute ExCluster statistics ####################################

    Stats_table <- ExClust_compute_stats(log2Clusters$gene, log2Clusters$pval)

    ######################################################################################################
    ###############################  FINAL DATA PREPARATION FOR OUTPUT ###################################
    ######################################################################################################

    ### ExCluster progression message
    message("Final ExCluster data organization & output ...","",sep="\n")

    # set up statistics & cluster columns in log2FC
    log2FC$pval <- NA
    log2FC$FDR <- NA
    log2FC$Cluster <- NA

    ### remove ENSG00000 from ENSEMBL IDs to make them easier to sort (for looping through each unique ID)
    # Now add an IDs column to the log2FC data frame
    log2FC$IDs <- parseGeneIDs(gene.Annot = rownames(log2FC))
    # now do the same for the IDs in the stats table
    Stats_table$IDs <- parseGeneIDs(Stats_table[,1])
    # now do the same for IDs in the log2Clusters
    log2Clusters$IDs <- parseGeneIDs(log2Clusters$gene)
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
    for (i in seq(2,nrow(log2FC))){

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
                log2FC$pval[seq(Start,Stop)] <- Stats_table$pvals[StatsCounter]
                log2FC$FDR[seq(Start,Stop)] <- Stats_table$FDRs[StatsCounter]
                StatsCounter <- min((StatsCounter+1), nrow(Stats_table))
            }

            ### Now check if the log2Clusters ID is equal to the current log2FC gene ID
            if (log2FC$IDs[Start] == log2Clusters$IDs[log2Start]){
                # if above is TRUE, grab data for the next 300 rows in log2Clusters
                # make sure we don't run out of rows in log2Start
                Upper <- min(300,(nrow(log2Clusters)-log2Start))
                TempData <- log2Clusters[seq(log2Start,(log2Start+Upper)),]
                # isolate data for only the current gene
                TempData <- TempData[which(TempData$IDs%in%TempData$IDs[1]),]
                # grab the exon bins from this gene
                ExonBins <- TempData$Bin
                # now loop through the exon bins in TempData & assign cluster numbers to log2FC
                for (j in seq(length(ExonBins))){
                    BinIndices <- which(log2FC$Bins[seq(Start,Stop)]%in%ExonBins[j])
                    log2FC$Cluster[seq(Start,Stop)][BinIndices] <- TempData$clustnum[j]
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
    final.log2FC <- data.frame(EnsG, ExonBins ,log2FC$Genes,log2FC$Cluster,log2FC$log2FC,
                               log2FC$log2Variance,log2FC$pval,log2FC$FDR,GFF.Chr,GFF.Strand,GFF.Start,GFF.End)
    # rename columns
    colnames(final.log2FC) <- c("EnsG","Exon_bin","Gene_name","Cluster",
                                "log2FC","log2Var","pval","FDR","Chr","Strand","Start","End")
    # re-add the read counts at the end of the dataframe
    final.log2FC <- data.frame(final.log2FC,log2FC[,seq(length(cond.Nums))])
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
            write.table(final.log2FC,file=paste(fullFilename,".txt",sep=""),sep="\t",
                        quote=FALSE,row.names=FALSE,col.names=TRUE)
        }else{
            warning(call = ExCluster_errors$ExClust_bad_outDir)
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
                dir.create(path = out.Dir,recursive = TRUE)
            }
            # make sure that FDR.cutoff is less than 0.2, or set it back to 0.05
            if (FDR.cutoff > 0.2){
                warning(call = ExCluster_errors$FDR_cutoff_too_high)
                FDR.cutoff <- 0.05
            }
            # run plotting function
            plotExonlog2FC(results.Data=final.log2FC, out.Dir=out.Dir, FDR.cutoff=FDR.cutoff)
        }else{
            warning(call = ExCluster_errors$plot_results_no_outDir)
        }
    }

    ### return final data frame resutls
    return(final.log2FC)
}
