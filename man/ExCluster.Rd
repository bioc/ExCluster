\name{ExCluster}
\alias{ExCluster}
\title{
ExCluster is the main function in this package, which clusters exon bins within each gene to determine signficiantly differentially expressed exon clusters.
}
\description{
The ExCluster function takes an input of normalized exon bin read counts per sample, and only accepts data from exactly two conditions at once, requiring two biological replicates per condition. ExCluster also requires GFF annotations to annotate exon bins
}
\usage{
ExCluster(exon.Counts, cond.Nums, annot.GFF, GFF.File, out.Dir,
    result.Filename, combine.Exons=TRUE, plot.Results=FALSE, FDR.cutoff=0.05)
}
\arguments{
  \item{exon.Counts}{exon.Counts must be assigned a data frame of normalized read counts per exon bin, which are attained from running the processCounts function. For example, if the results of processCounts are assigned to the normCounts variable, exon.Counts is assigned as: exon.Counts=normCounts This input is required.
    }
    \item{cond.Nums}{cond.Nums must be assigned an array of sample numbers corresponding, in exact order, to each BAM file counted and passed into the exon.Counts argument. The length of cond.Nums must also match the number of columns in exon.Counts exactly. For example, analyzing data with 3 biological replicates per each of two conditions, in that order, would be denoted as: cond.Nums=c(1,1,1,2,2,2) cond.Nums must be given at least 2 biological replicates per condition, and exactly 2 conditions -- otherwise, the ExCluster package will throw an error describing the problem. This is a required input.
    }
    \item{annot.GFF}{If this argument is specified, it must be given a GFF annotation data frame, which is generated from the GFF_convert function. For example, it can be specified as annot.GFF=GFF, assuming your GFF data frame is assigned to the GFF variable. Either annot.GFF or GFF.File must be specified. If both are specified, annot.GFF will take priority.
    }
    \item{GFF.File}{If this argument is specified, it must be given a full file path to a GFF file, including file name and extension. For example, this may look like: /Users/username/path/to/file.gff This argument is not required, but either annot.GFF or GFF.File must be specified.
    }
    \item{out.Dir}{The out.Dir argument should be assigned a character string specifying a full folder path where results may be written. This should be a specific folder to the data being analyzed, to avoid writing results from different analyses to the same folder.This argument is not required, but it is STRONGLY recommended that you write out the results of the ExCluster analysis, as this portion of the package can take well over 2 hours to complete. An example path for out.Dir would be: /Users/username/Documents/RNA-seq/project_A/
    }
    \item{result.Filename}{result.Filename should be given a character string specifying the specific filename for the ExCluster results table to be written to, within the out.Dir folder. File extension specification is not necessary, as '.txt.' will be added to the result.Filename regardless. This argument is not required, although it may be helpful to name your results by including both condition names, so the specific comparison made is easily identifiable at  a later point in time. By default, result.Filename will be assigned the value "ExClust_Results". An example filename that may be helpful could be: "ExClust_res_GeneA_shRNA_vs_scramble"
    }
    \item{combine.Exons}{combine.Exons should be assigned a logical value of either TRUE or FALSE. This denotes whether exon bins which are always co-expressed in the same transcripts should be combined into 'super-exons'. Doing this can be helpful to increase exon read depth, and it greatly reduces computation time. However, this should only be done in a standard RNA-seq analysis, when no instances of abberant splicing are predicted. If one suspects aberrant splicing in one of your conditions, this argument should not be set. By default combine.Exons=FALSE
    }
    \item{plot.Results}{plot.Results should be assigned a logical value of either TRUE or FALSE. This determines whether or not the ExCluster function should automatically run the plotExonlog2FC function, and plot exon bin log2FCs for each significantly differentially expressed gene. It is generally helpful to run this alongside your analysis, as it saves time. However, your ExCluster results can be saved and read back into R at a later date, from which you can run the plotExonlog2FC function separately. By default plot.Results is set to FALSE.
    }
    \item{FDR.cutoff}{The FDR.cutoff argument should be assigned a value between 0.01 and 0.2. Using FDR cutoffs outside these bounds is not recommended. This number determines which false discovery rate cutoff will be used to discover significant genes by ExCluster. However, this parameter is only used if plot.Results is specified. By default FDR.cutoff=0.05
    }
}
\examples{
# specify the path to the normCounts file in the ExCluster package
countsPath<- system.file("extdata","normCounts.txt",package="ExCluster")
# now read in the normCounts.txt file
normCounts <- read.table(file=countsPath,header=TRUE,row.names=1,
    stringsAsFactors=FALSE)
# now grab the path to the sub-sampled example GFF file
GFF_file <- system.file("extdata","sub_gen.v23.ExClust.gff3",package="ExCluster")
# assign condition numbers to your samples (we have 4 samples, 2 replicates per condition)
condNums <- c(1,1,2,2)
# now we run ExCluster, assigning its output to the ExClustResults variable
# we are not writing out the ExClustResults table, nor are we plotting exons
# we also use combine.Exons=TRUE, since we are conducting a standard analysis
ExClust_Results <- ExCluster(exon.Counts=normCounts,cond.Nums=condNums,
    GFF.File=GFF_file)
}
\value{
This is the main function of the ExCluster package, and returns a data frame of exon bin log2FCs, log2 read variances, exon clustering per gene, p-values, FDRs, and normalized exon counts. The results of the ExCluster function should typically be assigned to a variable, such as 'ExClustResults'.
}
\references{
Charrad M., Ghazzali N., Boiteau V., Niknafs A. (2014). NbClust: An R Package for Determining the Relevant Number of Clusters in a Data Set. Journal of Statistical Software, 61(6), 1-36.
}
\author{
R. Matthew Tanner
}
\note{
The ExCluster packages uses a scatter-distance index function to optimally cut hierarchically clustered exons within each gene. The code for the functions required were adapted in part, or in whole, from the NbClust R package (Charrad et al., 2014). These sections of the code are explicity specified, and the authors of NbClust provide no warranty for the use and fucntionality of said code.
}
