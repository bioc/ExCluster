GFF_convert <- function(GTF.File=NULL,GFF.File=NULL){

    ### Check to ensure that the GTF variable was given, GFF variable is optional
    if (is.null(GTF.File == TRUE)){
        stop(call. = "Error! You did not properly specify a full GTF filepath when running this function.
             A GTF file is required to proceed -- please supply a full GTF file path to the GTF.File argument.")
    }

    ### Now check to make sure that the GTF file path specified actually exists
    if (file.exists(GTF.File) == FALSE){
        stop(call. = "Error! The GTF file path that you specified does not exist. Please double check your GTF files full path & file name.
             For example, file paths should be entered as such: /User/username/path/to/file.gtf when given to the GTF.File argument.")
    }

    ### load GTF and set up final data frame
    temp.gtf.data<-read.table(GTF_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)

    ### now double check to make sure that 9 columns exactly are present in the GTF file
    if (ncol(temp.gtf.data) != 9){
        stop(call. = "Error! Your GTF file did not have 9 tab delimited columns. If you have edited your GTF file, please ensure the original number of columns and type of columns are retained.
             Additionally, if you have edited your GTF file, please ensure that tabs were used to delimit columns.
             If this error persists, please re-download your GTF annotations from the source database, and try again.")
    }

    # select only GTF rows that correspond to exons (includes UTRs, CDS, and retained introns)
    gtf.data<-temp.gtf.data[which(temp.gtf.data$V3 =='exon'),]
    # clean up temp.gtf.data
    rm(temp.gtf.data)


    # set up final plus and neg strand data frames for final output
    final.dataframe.plus.final<-as.data.frame(matrix(ncol=9,nrow=0))
    final.dataframe.neg.final<-as.data.frame(matrix(ncol=9,nrow=0))
    # assign plus and neg gtf.data to data frames
    gtf.plus<-gtf.data[which(gtf.data$V7== '+'),]
    gtf.neg<-gtf.data[which(gtf.data$V7== '-'),]
    # grab unique chromosome identifiers
    chrom.num<-unique(gtf.data$V1)
    # clean up gtf.data
    rm(gtf.data)

    if (nrow(gtf.plus) > 0){
        for (a in 1:length(chrom.num)){
            chrom.data <- gtf.plus[which(gtf.plus$V1 == chrom.num[a]),]
            final.dataframe.plus.chrm<-as.data.frame(matrix(ncol=9,nrow=0))

            genes.plus <- sub('.*gene_id ','',chrom.data$V9)
            genes.plus <- sub(';.*','',genes.plus)
            trans.plus <- sub('.*transcript_id ','',chrom.data$V9)
            trans.plus <- sub(';.*','',trans.plus)
            exons.org.plus <- sub('.*exon_number ','',chrom.data$V9)
            exons.org.plus <- paste("exon_number_",sub(';.*','',exons.org.plus),sep="")
            chrom.data$V10<-genes.plus
            chrom.data$V11<-trans.plus
            chrom.data$V12<-exons.org.plus
            genes=unique(genes.plus)

            ### index gene start & stops to make run faster
            Counter <- 1
            starts.plus <- array()
            stops.plus <- array()
            starts.plus[1] <- 1
            for (i in 2:length(genes.plus)){
                if (genes.plus[i] != genes.plus[(i-1)]){
                    stops.plus[Counter] <- i-1
                    Counter <- Counter + 1
                    starts.plus[Counter] <- i
                }
            }
            stops.plus[Counter] <- length(genes.plus)
            for (b in 1:length(genes)){
                #tempgene.data<-chrom.data[chrom.data$V10==b,]
                tempgene.data <- chrom.data[starts.plus[b]:stops.plus[b],]

                ### grab corresponding gene names
                GeneNames <- sub('.*gene_name ','',tempgene.data$V9)
                GeneNames <- sub(';.*','',GeneNames)
                GeneName <- GeneNames[1]

                exon.range<-c()
                ### generate exon and transcript ranges for ID of transcripts per exon bin
                transcript.ranges<-list()
                for (c in 1:nrow(tempgene.data)){
                    exon.range<-c(exon.range,tempgene.data[c,4]:tempgene.data[c,5])
                    transcript.ranges[[length(transcript.ranges)+1]]<-c(tempgene.data[c,4]:tempgene.data[c,5])
                }
                exon.range<-sort(unique(exon.range))
                diff.range<-diff(exon.range)
                diff.pos<-which(diff.range>1)
                gaps<-c(max(exon.range),exon.range[diff.pos])
                exon.edges<-sort(unique((c(tempgene.data[,4],tempgene.data[,5]))))
                Counter <- 1
                exon.start <- array()
                exon.end <- array()
                #exon.start<-c(gene.min)
                #exon.end<-c(gene.max)

                for (d in 1:length(exon.edges)){
                    if (is.na(match(exon.edges[d],gaps))== TRUE){
                        #exon.start<-c(exon.start,exon.edges[d])
                        #exon.end<-c(exon.end,exon.edges[d+1]-1)
                        exon.start[Counter] <- exon.edges[d]
                        exon.end[Counter] <- (exon.edges[d+1]-1)
                        Counter <- Counter + 1
                    }else{
                        #exon.end[length(exon.end)]<-exon.end[length(exon.end)]+1
                        exon.end[Counter-1] <- (exon.end[length(exon.end)]+1)
                    }
                }
                gene.dataframe<-as.data.frame(matrix(ncol=9,nrow=length(exon.start)))
                chr.collumn<-rep(chrom.num[a],length(exon.start))
                bin.nums <- c(1:length(exon.start))
                genes.collumn<-sprintf(paste(genes[b],":%03d",sep=""),bin.nums)
                type.collumn<-c(paste("exonic_part_",1:(length(exon.start)),sep=""))
                dot.collumn<-rep(".",length(exon.start))
                plus.collumn<-rep("+",length(exon.start))
                names.column <- rep(GeneName,length(exon.start))

                ### generate master whole exon ranges
                master.ranges <- list()
                for (i in 1:nrow(tempgene.data)){
                    master.ranges[[i]] <- tempgene.data$V4[i]:tempgene.data$V5[i]
                }

                ### now loop through each exon bin 'exon.start' and through each master exon to find overlap
                ### transcripts to assign to each exon bin
                bin.trans <- array()
                for(e in 1:length(exon.start)){
                    Counter <- 1
                    temp.trans <- array()
                    for (m in 1:length(master.ranges)){
                        if ((exon.start[e]%in%master.ranges[[m]]) == TRUE){
                            temp.trans[Counter] <- tempgene.data$V11[m]
                            Counter <- Counter+1
                        }
                    }
                    bin.trans[e] <- paste(unique(temp.trans),collapse="+")
                }

                gene.dataframe$V1<-chr.collumn
                gene.dataframe$V2<-genes.collumn
                gene.dataframe$V3<-type.collumn
                gene.dataframe$V4<-exon.start
                gene.dataframe$V5<-exon.end
                gene.dataframe$V6<-dot.collumn
                gene.dataframe$V7<-plus.collumn
                gene.dataframe$V8<-dot.collumn
                gene.dataframe$V9<-bin.trans
                gene.dataframe$V10<-names.column

                final.dataframe.plus.chrm<-rbind(final.dataframe.plus.chrm,gene.dataframe)
            }
            position.dataframe.plus<-cbind(final.dataframe.plus.chrm$V4,final.dataframe.plus.chrm$V5)
            duplicate.tracker<-which(duplicated(position.dataframe.plus) | duplicated(position.dataframe.plus, fromLast = TRUE) == TRUE)
            if (length(duplicate.tracker) > 0 ){
                final.dataframe.plus.final<-rbind(final.dataframe.plus.final, final.dataframe.plus.chrm[-duplicate.tracker,])
            }else{
                final.dataframe.plus.final<-rbind(final.dataframe.plus.final, final.dataframe.plus.chrm)
            }
        }
    }

    if (nrow(gtf.neg) > 0){
        for (a in 1:length(chrom.num)){
            chrom.data <- gtf.neg[which(gtf.neg$V1 == chrom.num[a]),]
            final.dataframe.neg.chrm<-as.data.frame(matrix(ncol=9,nrow=0))

            ### break "meta data" column 9 into columns, sep = ;
            #out <- strsplit(as.character(gtf.neg$V9),';')
            #metaCol <- do.call(rbind, out)
            genes.neg <- sub('.*gene_id ','',chrom.data$V9)
            genes.neg <- sub(';.*','',genes.neg)
            trans.neg <- sub('.*transcript_id ','',chrom.data$V9)
            trans.neg <- sub(';.*','',trans.neg)
            exons.org.neg <- sub('.*exon_number ','',chrom.data$V9)
            exons.org.neg <- paste("exon_number_",sub(';.*','',exons.org.neg),sep="")
            chrom.data$V10<-genes.neg
            chrom.data$V11<-trans.neg
            chrom.data$V12<-exons.org.neg
            genes=unique(genes.neg)

            ### index gene start & stops to make run faster
            Counter <- 1
            starts.neg <- array()
            stops.neg <- array()
            starts.neg[1] <- 1
            for (i in 2:length(genes.neg)){
                if (genes.neg[i] != genes.neg[(i-1)]){
                    stops.neg[Counter] <- i-1
                    Counter <- Counter + 1
                    starts.neg[Counter] <- i
                }
            }
            stops.neg[Counter] <- length(genes.neg)


            for (b in 1:length(genes)){
                tempgene.data <- chrom.data[starts.neg[b]:stops.neg[b],]

                ### grab corresponding gene names
                GeneNames <- sub('.*gene_name ','',tempgene.data$V9)
                GeneNames <- sub(';.*','',GeneNames)
                GeneName <- GeneNames[1]

                exon.range<-c()
                ### generate exon and transcript ranges for ID of transcripts per exon bin
                transcript.ranges<-list()

                for (c in 1:nrow(tempgene.data)){
                    exon.range<-c(exon.range,tempgene.data[c,4]:tempgene.data[c,5])
                    transcript.ranges[[length(transcript.ranges)+1]]<-c(tempgene.data[c,4]:tempgene.data[c,5])
                }
                exon.range<-sort(unique(exon.range))
                diff.range<-diff(exon.range)
                diff.pos<-which(diff.range>1)
                gaps<-c(max(exon.range),exon.range[diff.pos])
                exon.edges<-sort(unique((c(tempgene.data[,4],tempgene.data[,5]))))
                Counter <- 1
                exon.start <- array()
                exon.end <- array()

                #exon.start<-c(gene.min)
                #exon.end<-c(gene.max)

                for (d in 1:length(exon.edges)){
                    if (is.na(match(exon.edges[d],gaps))== TRUE){
                        #exon.start<-c(exon.start,exon.edges[d])
                        #exon.end<-c(exon.end,exon.edges[d+1]-1)
                        exon.start[Counter] <- exon.edges[d]
                        exon.end[Counter] <- (exon.edges[d+1]-1)
                        Counter <- Counter + 1
                    }else{
                        #exon.end[length(exon.end)]<-exon.end[length(exon.end)]+1
                        exon.end[Counter-1] <- (exon.end[length(exon.end)]+1)
                    }
                }
                gene.dataframe<-as.data.frame(matrix(ncol=9,nrow=length(exon.start)))
                chr.collumn<-rep(chrom.num[a],length(exon.start))
                bin.nums <- c(length(exon.start):1)
                genes.collumn<-sprintf(paste(genes[b],":%03d",sep=""),bin.nums)
                type.collumn<-c(paste("exonic_part_",(length(exon.start):1),sep=""))
                dot.collumn<-rep(".",length(exon.start))
                neg.collumn<-rep("-",length(exon.start))
                names.column <- rep(GeneName,length(exon.start))

                ### generate master whole exon ranges
                master.ranges <- list()
                for (i in 1:nrow(tempgene.data)){
                    master.ranges[[i]] <- tempgene.data$V4[i]:tempgene.data$V5[i]
                }

                ### now loop through each exon bin 'exon.start' and through each master exon to find overlap
                ### transcripts to assign to each exon bin
                bin.trans <- array()
                for(e in 1:length(exon.start)){
                    Counter <- 1
                    temp.trans <- array()
                    for (m in 1:length(master.ranges)){
                        if ((exon.start[e]%in%master.ranges[[m]]) == TRUE){
                            temp.trans[Counter] <- tempgene.data$V11[m]
                            Counter <- Counter+1
                        }
                    }
                    bin.trans[e] <- paste(unique(temp.trans),collapse="+")
                }

                gene.dataframe$V1<-chr.collumn
                gene.dataframe$V2<-genes.collumn
                gene.dataframe$V3<-type.collumn
                gene.dataframe$V4<-exon.start
                gene.dataframe$V5<-exon.end
                gene.dataframe$V6<-dot.collumn
                gene.dataframe$V7<-neg.collumn
                gene.dataframe$V8<-dot.collumn
                gene.dataframe$V9<-bin.trans
                gene.dataframe$V10 <- names.column

                final.dataframe.neg.chrm<-rbind(final.dataframe.neg.chrm,gene.dataframe)
            }
            position.dataframe.neg<-cbind(final.dataframe.neg.chrm$V4,final.dataframe.neg.chrm$V5)
            duplicate.tracker<-which(duplicated(position.dataframe.neg) | duplicated(position.dataframe.neg, fromLast = TRUE) == TRUE)
            if (length(duplicate.tracker) > 0 ){
                final.dataframe.neg.final<-rbind(final.dataframe.neg.final, final.dataframe.neg.chrm[-duplicate.tracker,])
            }else{
                final.dataframe.neg.final<-rbind(final.dataframe.neg.final, final.dataframe.neg.chrm)
            }
        }
    }
    ### now combine positive and negative data
    final.dataframe<-rbind(final.dataframe.plus.final,final.dataframe.neg.final)
    # clean up just incase
    rm(final.dataframe.plus.final)
    rm(final.dataframe.neg.final)

    ### if the GFF outpath was specified, AND the filepath is writeable, write out the file
    # check if file path is writeable

    if (is.null(GFF.File) == FALSE){
        WriteCheck <- file.access(dirname(GFF.File), mode=2)
        if (WriteCheck == 0){
            write.table(final.dataframe, file=GFF.File,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
        }
    }

    return(final.dataframe)
}
