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
            for (n in 1:NRows){
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


    IndexGeneStartStop <- function(x){
        Counter <- 1
        gene.names <- NULL
        gene.starts <- NULL
        gene.stops <- NULL
        gene.starts[1] <- 1
        MAXrow <- length(x)
        for (i in 1:length(x)){
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

    #########################################################################################################
    ########################################### MAIN CODE ###################################################
    #########################################################################################################

    ### load GTF and set up final data frame
    gtf.data<-read.table(GTF.File, header=FALSE, sep="\t", stringsAsFactors=FALSE)

    ### now double check to make sure that 9 columns exactly are present in the GTF file
    if (ncol(gtf.data) != 9){
        stop(call = "Your GTF file did not have 9 tab delimited columns.
If you have edited your GTF file, please ensure the original number of columns and type of columns are retained.
Additionally, if you have edited your GTF file, please ensure that tabs were used to delimit columns.
If this error persists, please re-download your GTF annotations from the source database, and try again.")
    }

    # select only GTF rows that correspond to exons (includes UTRs, CDS, and retained introns)
    gtf.data<-gtf.data[which(gtf.data$V3 =='exon'),]
    # add columns to gtf.file for gene_name, gene_id, and transcript_id
    gtf.data$gene.id <- gsub(".*gene_id (.*?);.*", "\\1", gtf.data$V9)
    gtf.data$gene.name <- gsub(".*gene_name (.*?);.*", "\\1", gtf.data$V9)
    gtf.data$tx.id <- gsub(".*transcript_id (.*?);.*", "\\1", gtf.data$V9)
    # now remove that mess of a 9th column
    gtf.data <- gtf.data[,-c(9)]

    # assign plus and neg gtf.data to data frames
    gtf.plus<-gtf.data[which(gtf.data$V7== '+'),]
    gtf.neg<-gtf.data[which(gtf.data$V7== '-'),]
    rm(gtf.data)

    # make empty data frames for GFF output, in case one strand is not present
    plus.gff <-as.data.frame(matrix(ncol=10,nrow=0))
    neg.gff <-as.data.frame(matrix(ncol=10,nrow=0))


    ###### Process +ve strand genes first
    if (nrow(gtf.plus) > 0){
        RStart <- 1
        plus.gff <- matrix(NA,nrow=nrow(gtf.plus)*10,ncol=6)
        # generate a hash table for the start/stop of each gene (2nd and 3rd columns, 1st column gene_id)
        GeneHash <- IndexGeneStartStop(gtf.plus$gene.id)
        # now loop through each gene in the GeneHash table & collapse exon bins with an lapply function

        for (x in 1:nrow(GeneHash)){
            # subset the current gene data
            temp.data <- gtf.plus[c(GeneHash[x,2]:GeneHash[x,3]),]
            # collapse exon bins for this gene
            gff.data <- CollapseExons(temp.data)
            # exon bin numbers for the current gene (1:n for + strand)
            Bins <- 1:nrow(gff.data)
            # now add in exon bin numbers based on Bins
            gff.data[,2]<- sprintf(paste(gff.data[,2],":%03d",sep=""),Bins)
            # now add the data to the plus.gff matrix
            RStop <- RStart+length(Bins)-1
            plus.gff[RStart:RStop,] <- gff.data
            RStart <- RStop + 1
        }
        # now flatten the plus.gff list into a data.frame
        plus.gff <- plus.gff[1:RStop,]
        # add in the "exonic_part" column
        GFF.part <- paste("exonic_part_",sub(".*:","",plus.gff[,2]),sep="")
        # add in columns that are generic across all chromosomes in the + strand
        GFF.dot1 <- rep(".",nrow(plus.gff))
        GFF.strand <- rep("+",nrow(plus.gff))
        GFF.dot2 <- GFF.dot1
        plus.gff <- data.frame(plus.gff[,c(1:2)],GFF.part,plus.gff[,c(3:4)],GFF.dot1,GFF.strand,GFF.dot2,plus.gff[,c(5:6)])
        # remove the original gtf
        rm(gtf.plus)
    }

    ###### Now process -ve strand
    if (nrow(gtf.neg) > 0){
        RStart <- 1
        neg.gff <- matrix(NA,nrow=nrow(gtf.neg)*10,ncol=6)
        # generate a hash table for the start/stop of each gene (2nd and 3rd columns, 1st column gene_id)
        GeneHash <- IndexGeneStartStop(gtf.neg$gene.id)
        # now loop through each gene in the GeneHash table & collapse exon bins with an lapply function

        for (x in 1:nrow(GeneHash)){
            # subset the current gene data
            temp.data <- gtf.neg[c(GeneHash[x,2]:GeneHash[x,3]),]
            # collapse exon bins for this gene
            gff.data <- CollapseExons(temp.data)
            # exon bin numbers for the current gene (1:n for + strand)
            Bins <- nrow(gff.data):1
            # now add in exon bin numbers based on Bins
            gff.data[,2]<- sprintf(paste(gff.data[,2],":%03d",sep=""),Bins)
            # now add the data to the neg.gff matrix
            RStop <- RStart+length(Bins)-1
            neg.gff[RStart:RStop,] <- gff.data
            RStart <- RStop + 1
        }
        # now flatten the neg.gff list into a data.frame
        neg.gff <- neg.gff[1:RStop,]
        # add in the "exonic_part" column
        GFF.part <- paste("exonic_part_",sub(".*:","",neg.gff[,2]),sep="")
        # add in columns that are generic across all chromosomes in the - strand
        GFF.dot1 <- rep(".",nrow(neg.gff))
        GFF.strand <- rep("-",nrow(neg.gff))
        GFF.dot2 <- GFF.dot1
        neg.gff <- data.frame(neg.gff[,c(1:2)],GFF.part,neg.gff[,c(3:4)],GFF.dot1,GFF.strand,GFF.dot2,neg.gff[,c(5:6)])
        # remove the original gtf
        rm(gtf.neg)
    }

    # combine +ve and -ve strand data
    final.gff <- rbind(plus.gff,neg.gff)
    # name columns for sorting
    colnames(final.gff) <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10")
    # clean up
    rm(plus.gff)
    rm(neg.gff)
    # sort on position and then chromsome
    final.gff <- final.gff[order(final.gff[,1],as.numeric(as.character(final.gff[,4])),final.gff[,7]),]

    ### if the GFF outpath was specified, AND the filepath is writeable, write out the file
    # check if file path is writeable

    if (is.null(GFF.File) == FALSE){
        WriteCheck <- file.access(dirname(GFF.File), mode=2)
        if (WriteCheck == 0){
            write.table(final.gff, file=GFF.File,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
        }
    }
    return(final.gff)
}
