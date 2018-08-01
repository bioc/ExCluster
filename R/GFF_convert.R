GFF_convert <- function(GTF.File=NULL,GFF.File=NULL){

    ### make sure R doesn't add factors -- R creators never should have made this default = TRUE
    options(stringsAsFactors=FALSE)

    ### Check to ensure that the GTF variable was given, GFF variable is optional
    if (is.null(GTF.File == TRUE)){
        stop(call = "You did not properly specify a full GTF filepath when running this function.
    A GTF file is required to proceed -- please supply a full GTF file path to the GTF.File argument.")
    }

    ### Now check to make sure that the GTF file path specified actually exists
    if (file.exists(GTF.File) == FALSE){
        stop(call = "The GTF file path that you specified does not exist.
    Please double check your GTF file's full path & file name.
    For example, a GTF.File path for mac/linux users might look like:
        /User/username/path/to/file.gtf")
    }

    #########################################################################################################
    ########################################### FUNCTIONS ###################################################
    #########################################################################################################


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

    ################################# Check GTF file for correct formatting #################################

    ### load GTF and run some checks to verify the file format integrity
    gtf.data<-read.table(GTF.File, header=FALSE, sep="\t", stringsAsFactors=FALSE)

    ### double check to make sure that 9 columns exactly are present in the GTF file
    if (ncol(gtf.data) != 9){
        stop(call = "Your GTF file did not have 9 tab delimited columns.
If you have edited your GTF file, please ensure the original number of columns and type of columns are retained.
Additionally, if you have edited your GTF file, please ensure that tabs were used to delimit columns.
If this error persists, please re-download your GTF annotations from the source database, and try again.")
    }

    ### check to make sure column 3 has values 'gene' 'exon' and 'transcript'
    feature.Values <- unique(gtf.data[,3])
    # check to see if "gene", "transcript", and "exon" are present in feature.Values
    feature.Indices <- match(c("gene","transcript","exon"),feature.Values)
    # the length of the above should be 3
    if (length(feature.Indices) != 3){
        stop(call = "Your GTF file did not have 'gene', 'transcript', and 'exon' features present in column 3.
If you have edited your GTF file, please ensure the original number of columns and type of columns are retained.
Also double check that your GTF file was read into R without row or column names, as the GTF format has none.
If this error persists, please re-download your GTF annotations from the source database, and try again.")
    }

    ### check to make sure columns 4 & 5 are numeric start/stops with positive differences only
    # we do this by subtracting 'start' from 'stop' and making sure the math works, and minimum = 0
    # test to see if the math will work (will only work if columns 4 & 5 are properly formatted)
    coord.Diff <- try(as.numeric(gtf.data[,5]) - as.numeric(gtf.data[,4]),silent = TRUE)
    # check to see if we have an error
    if (substr(coord.Diff[1],1,5) == "Error"){
        stop(call = "The start and stop columns (4 & 5) of your GTF file were not properly formatted.
Please ensure your gtf.file was read correctly, without row or column names from the original GTF file.
If this error persists, please re-download your GTF annotations from the source database, and try again.")
    }
    # if that did not fail, the code continues and we verify that all stop - start subtractions are >= 0
    if (min(coord.Diff) < 0){
        stop(call = "There were negative distances between your stop and start locations.
The start/stop locations of your genomic features in this GTF file appear to be corrupted.
Please ensure column 4 of your GTF file is start location, and column 5 is stop location.
If this error persists, please re-download your GTF annotations from the source database, and try again.")
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
        stop(call = "The strand column of your GTF file (column 7) had more than 2 unique values.
This column should only have '+' and '-' values. Please ensure your GTF file is read into R without row/column names.
If this error persists, please re-download your GTF annotations from the source database, and try again.")
    }
    # if the previous check passed, make sure at least one of '-' or '+' is present
    if ((positive.Length + negative.Length) < 1){
        stop(call = "The strand column of your GTF file (column 7) did not contain '+' or '-' values.
This column should only have '+' and/or '-' values. Please ensure your GTF file is read into R without row/column names.
If this error persists, please re-download your GTF annotations from the source database, and try again.")
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
        stop(call = "Column 9 in your GTF file did not contain gene_id and/or transcript_id values in the correct places.
Please verify that your GTF file was read into R without row or column names, and that each row has a gene_id in column 9.
If this error persists, please re-download your GTF annotations from the source database, and try again.")
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
        temp.gff <- matrix(NA,nrow=nrow(exons.collapsed),ncol=10)
        # add chromosome names
        temp.gff[,1] <- exons.collapsed[,1]
        # add EnsID + exon bin
        temp.gff[,2] <- sprintf(paste(exons.collapsed[,2],":%03d",sep=""),Bins)
        # add exon_bin column
        temp.gff[,3] <- sprintf("exonic_part_%03d",Bins)
        # add start + stop columns
        temp.gff[,c(4,5)] <- exons.collapsed[,c(3,4)]
        # add '.' columns
        temp.gff[,c(6,8)] <- "."
        # add strand column (+ve or -ve based on previous chr.Strand)
        temp.gff[,7] <- chr.Strand
        # add transcript column & gene name column
        temp.gff[,c(9,10)] <- exons.collapsed[,c(5,6)]
        # assign temp.gff to be item x on the gff.data list
        gff.data[[x]] <- temp.gff
    }
    # clean up gtf data now
    rm(gtf.data)
    # now we can combine the gff data into one data frame structure with do.call
    gff.data <- data.frame(do.call(rbind,gff.data))
    # change filenames to arbitrary V1, V2, etc. (necessary for next step in read counting)
    colnames(gff.data) <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10")
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
    return(gff.data)
}
