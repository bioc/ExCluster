\name{plotExonlog2FC}
\alias{plotExonlog2FC}
\title{
The plotExonlog2FC function plots exon log2FC means and variances for significantly differentially expressed genes in ExCluster results.
}
\description{
This function takes an input of ExCluster's results dataframe, as well as a directory to write files to, and plots exon bin log2FC means and variances for each significantly differentially expressed gene. Significant genes are determined by applying an FDR cutoff to the ExCluster results. By default, FDR.cutoff=0.05, although this parameter may be set between 0.01 and 0.2.
}
\usage{
plotExonlog2FC(results.Data, out.Dir, FDR.cutoff)
}
\arguments{
  \item{results.Data}{This argument should be assigned a data frame containing the results from an ExCluster analysis, obtained by running the function ExCluster. This argument is required to be set. Assuming ExCluster results are assigned to the variable clustResults, this parameter woud be set as follows: results.Data=clustResults
    }
    \item{out.Dir}{The outDir argument should be assigned a character string specifying a full folder path where exon bin images per gene may be written. To keep results coherent, this should be a sub-folder within the folder the original ExCluster results were written to. This directory need not exist, as R will recursively create it. However, this argument is required. For example, if ExCluster results were written to the folder /Users/username/Documents/RNA-seq/project_A/ then the folder the images should be something like: /Users/username/Documents/RNA-seq/project_A/exon_log2FC_images/
    }
    \item{FDR.cutoff}{The FDR.cutoff argument should be assigned a value between 0.01 and 0.2. Using FDR cutoffs outside these bounds is not recommended. This number determines which false discovery rate cutoff will be used to discover significant genes by ExCluster. This parameter is not required, and by default FDR.cutoff=0.05
    }
}
\examples{
# grab the file path to sub-sampled results from an ExCluster function analysis
ExClust_respath <- system.file("extdata","sub_ExClust_results.txt",package="ExCluster")
# load this data into R
ExClust_Results <- read.table(file=ExClust_respath,sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
# now we must specify a directory to write images to
outDir <- paste(system.file(package="ExCluster"),"/extdata/Images/",sep="")
# now we can run our exon log2FC plotting function
plotExonlog2FC(results.Data=ExClust_Results, out.Dir=outDir)
}
\value{
This function does not return any values, but instead plots exon log2FC means and variances for significant results from an ExCluster function results data.frame.
}
\author{
R. Matthew Tanner
}
