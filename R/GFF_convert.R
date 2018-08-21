GFF_convert <- function(annot.GTF=NULL, GTF.File=NULL,GFF.File=NULL){

    ### make sure R doesn't add factors -- R creators never should have made this default = TRUE
    options(stringsAsFactors=FALSE)

    ### Check to ensure that the GTF variable was given, GFF variable is optional
    if (is.null(GTF.File) == TRUE){
        if (is.null(annot.GTF) == TRUE){
            stop(call = ExCluster_errors$GTF.File_missing)
        }
    }

    ### Now check that the GTF file path specified actually exists if annot.GTF == false
    if (is.null(annot.GTF) == TRUE){
        if (file.exists(GTF.File) == FALSE){
            stop(call = ExCluster_errors$bad_GTF_filepath)
        }else{
            # if the GTF.File exists, read it into R
            gtf.data<-read.table(GTF.File, header=FALSE, sep="\t", stringsAsFactors=FALSE)
        }
    }else{

        ### if annot.GTF exists, read it into R
        # if the typeof(annot.GTF) is S4, treat it as an rtracklayer object
        if (typeof(annot.GTF) == "S4"){
            # reformat rtracklayer to GTF
            gtf.data <- reformat_GTF(rtracklayer.GTF=annot.GTF)
        }else{
            # rename annot.GTF
            gtf.data <- annot.GTF
        }
        # clean up
        rm(annot.GTF)
    }

    ######################## Check GTF file for correct formatting ########################

    ### run some checks to verify the GTF file format integrity

    ### double check to make sure that 9 columns exactly are present in the GTF file
    if (ncol(gtf.data) != 9){
        stop(call = ExCluster_errors$check_GTF_columns)
    }

    ### check to make sure column 3 has values 'gene' 'exon' and 'transcript'
    feature.Values <- unique(gtf.data[,3])
    # check to see if "gene", "transcript", and "exon" are present in feature.Values
    feature.Indices <- match(c("gene","transcript","exon"),feature.Values)
    # the length of the above should be 3
    if (length(feature.Indices) != 3){
        stop(call = ExCluster_errors$GTF_missing_features)
    }

    ### check to make sure columns 4 & 5 are numeric start/stops with positive differences only
    # we do this by subtracting 'start' from 'stop' and making sure the math works, and minimum = 0
    # test to see if the math will work (will only work if columns 4 & 5 are properly formatted)
    coord.Diff <- try(as.numeric(gtf.data[,5]) - as.numeric(gtf.data[,4]),silent = TRUE)
    # check to see if we have an error
    if (substr(coord.Diff[1],1,5) == "Error"){
        stop(call = ExCluster_errors$bad_GTF_start_stop)
    }
    # if that did not fail, the code continues and we verify that all stop - start subtractions are >= 0
    if (min(coord.Diff) < 0){
        stop(call = ExCluster_errors$start_stop_not_numeric)
    }

    ### check to make sure column 7 unique values have at least "+" or "-" (could be just one for test gtf files)
    # unique values of column 7 in the GTF file
    strand.Values <- unique(gtf.data[,7])
    # is "+" present? if true this will == 1, if false this will == 0
    positive.Length <- length(which('+'%in%strand.Values))
    # is "-" present? if true this will == 1, if false this will == 0
    negative.Length <- length(which('-'%in%strand.Values))
    # make sure there are no more than 2 unique values of column 7
    if (length(strand.Values) > 2){
        stop(call = ExCluster_errors$bad_GTF_plus_minus)
    }
    # if the previous check passed, make sure at least one of '-' or '+' is present
    if ((positive.Length + negative.Length) < 1){
        stop(call = ExCluster_errors$GTF_plus_minus_missing)
    }

    ### lastly, check to make sure transcript_id and gene_id values are present in GTF metadata column 9
    # grab 1st 25 rows -- sufficient for testing purposes
    metacol.Values <- gtf.data[1:25,9]
    # check the number of rows (if any) 'gene_id' is present in for this metadata column
    gene_id.Length <- length(grep("gene_id",metacol.Values))
    # check the number of rows (if any) 'transcript_id' is present in for this metadata column
    transcript_id.Length <- length(grep("transcript_id",metacol.Values))
    # now make sure gene_id.Length == 25 and transcript_id.Length > 1
    # this is because each row should have a gene_id, and at least some rows should contain a transcript_id
    if (gene_id.Length != 25 || transcript_id.Length < 1){
        stop(call = ExCluster_errors$GTF_missing_gene_id)
    }

    #########################################################################################################
    ########################################### MAIN CODE ###################################################
    #########################################################################################################

    # select only GTF rows that correspond to exons (includes UTRs, CDS, and retained introns)
    gtf.data<-gtf.data[which(gtf.data$V3 =='exon'),]
    # add columns to gtf.file for gene_name, gene_id, and transcript_id
    gtf.data$gene.id <- gsub(".*gene_id (.*?);.*", "\\1", gtf.data$V9)
    gtf.data$gene.name <- gsub(".*gene_name (.*?);.*", "\\1", gtf.data$V9)
    gtf.data$tx.id <- gsub(".*transcript_id (.*?);.*", "\\1", gtf.data$V9)
    # now remove that mess of a 9th column
    gtf.data <- gtf.data[,-c(9)]
    # make an empty data matrix for gff file
    gff.data<- vector("list", 1)
    # generate a hash table for the start/stop of each gene (2nd and 3rd columns, 1st column gene_id)
    GeneHash <- IndexGeneStartStop(gtf.data$gene.id)
    # now loop through each gene in the GeneHash table & collapse exon bins with an lapply function
    for (x in seq(nrow(GeneHash))){
        # subset the current gene data
        temp.data <- gtf.data[c(GeneHash[x,2]:GeneHash[x,3]),]
        # collapse exon bins for this gene
        exons.collapsed <- CollapseExons(temp.data)
        # exon bin numbers for the current gene (1:n for +ve strand by default, n:1 for -ve strand)
        Bins <- seq(nrow(exons.collapsed))
        # strand is by default +ve
        chr.Strand <- "+"
        # if -ve strand, reverse this
        if (temp.data[1,7] == "-"){
            Bins <- rev(Bins)
            # also change chr.Strand to -ve
            chr.Strand <- "-"
        }
        # make temporary gff matrix
        temp.gff <- matrix(NA,nrow=nrow(exons.collapsed),ncol=9)
        # add chromosome names
        temp.gff[,1] <- exons.collapsed[,1]
        # add EnsID + exon bin
        temp.gff[,2] <- "."
        # add exon_bin column
        temp.gff[,3] <- "exon"
        # add start + stop columns
        temp.gff[,c(4,5)] <- exons.collapsed[,c(3,4)]
        # add '.' columns
        temp.gff[,c(6,8)] <- "."
        # add strand column (+ve or -ve based on previous chr.Strand)
        temp.gff[,7] <- chr.Strand
        # prepare gene IDs, Names, and Transcripts
        gff.ID <- sprintf(paste("ID=",exons.collapsed[,2],":%03d",sep=""),Bins)
        gff.Name <- sprintf(paste("Name=",exons.collapsed[,6],sep=""))
        gff.Transcript <- paste("Transcripts=",exons.collapsed[,5],sep="")
        # add a final 9th column with gene IDs, names, and transcripts
        temp.gff[,9] <- paste(gff.ID,gff.Name,gff.Transcript,sep=";")
        rm(gff.ID)
        rm(gff.Name)
        rm(gff.Transcript)
        # assign temp.gff to be item x on the gff.data list
        gff.data[[x]] <- temp.gff
    }
    # clean up gtf data now
    rm(gtf.data)
    # now we can combine the gff data into one data frame structure with do.call
    gff.data <- data.frame(do.call(rbind,gff.data))
    # change filenames to arbitrary V1, V2, etc. (necessary for next step in read counting)
    colnames(gff.data) <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9")
    # sort on position and then chromsome
    gff.data <- gff.data[order(gff.data[,1],as.numeric(as.character(gff.data[,4])),gff.data[,7]),]

    ### if the GFF outpath was specified, AND the filepath is writeable, write out the file
    # check if file path is writeable
    if (is.null(GFF.File) == FALSE){
        WriteCheck <- file.access(dirname(GFF.File), mode=2)
        if (WriteCheck == 0){
            write.table(gff.data, file=GFF.File,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
        }
    }

    # reformat GFF3 data to GRanges object
    gff.data <- GRangesFromGFF(gff.data)
    # return & end function
    return(gff.data)
}
