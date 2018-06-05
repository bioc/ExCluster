### This R script contains the ExCluster vignette R executable code for testing ExCluster

###############################################################
### Code chunk #1 for ExCluster.Rnw vignette: lines 183-188 ###
###############################################################

library(ExCluster)
# load the sub-sampled GTF file path from the ExCluster package
GTF_file <- system.file("extdata","sub_gen.v23.gtf", package = "ExCluster")
# now run GTF_file without assigning a GFF_file to write out, assigning the results to the GFF object
GFF <- GFF_convert(GTF.File=GTF_file)

###############################################################
### Code chunk #2 for ExCluster.Rnw vignette: lines 195-206 ###
###############################################################

library(ExCluster)
# specify the path to the ExCluster package
ExClust_Path <- system.file(package="ExCluster")
# now find the bam files within that folder
bamFiles <- list.files(ExClust_Path,recursive=TRUE,pattern="*.bam",full.names=TRUE)
# now grab the path to the sub-sampled example GFF file
GFF_file <- system.file("extdata","sub_gen.v23.ExClust.gff",package="ExCluster")
# assign sample names (only 2 replicates per condition in this example)
sampleNames <- c("iPSC_cond1_rep1","iPSC_cond1_rep2","iPSC_cond2_rep1","iPSC_cond2_rep2")
# now run processCounts, with pairedReads=TRUE because we are counting paired-end data
normCounts <- processCounts(bam.Files=bamFiles, sample.Names=sampleNames, GFF.File=GFF_file, pairedReads=TRUE)

###############################################################
### Code chunk #3 for ExCluster.Rnw vignette: lines 233-242 ###
###############################################################

library(ExCluster)
# specify the path to the normCounts file in the ExCluster package
countsPath<- system.file("extdata","normCounts.txt",package="ExCluster")
# now read in the normCounts.txt file
normCounts <- read.table(file=countsPath,header=TRUE,row.names=1,stringsAsFactors=FALSE)
# now grab the path to the sub-sampled example GFF file
GFF_file <- system.file("extdata","sub_gen.v23.ExClust.gff",package="ExCluster")
# assign condition numbers to your samples (we have 4 samples, 2 replicates per condition)
condNums <- c(1,1,2,2)
# now we run ExCluster, assigning its output to the ExClustResults variable
# we are not writing out the ExClustResults table, nor are we plotting exons
# we also use CombineExons=TRUE, since we are conducting a standard analysis
ExClust_Results <- ExCluster(exonCounts=normCounts,cond.Nums=condNums,GFF.File=GFF_file,CombineExons=TRUE)


