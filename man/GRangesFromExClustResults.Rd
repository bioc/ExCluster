\name{GRangesFromExClustResults}
\alias{GRangesFromExClustResults}
\title{
Converts ExCluster results data table to GRanges format.
}
\description{
GRangesFromExClustResults converts the results from the main ExCluster function into GRanges objects using the GRanges package.
}
\usage{
GRangesFromExClustResults(results.ExClust)
}
\arguments{
    \item{results.ExClust}{This argument should be given the table of data that the main ExCluster function produces. This table of data should contain at minimum 12 columns, plus exin bin read count columns if they are present.
    }
}
\examples{
# grab the path to the sub-sampled ExCluster toy dataset results from the ExCluster function
ExClust_Results <- system.file("extdata","sub_ExClust_results.txt",package="ExCluster")
# read this data into R
ExClust_Results <- read.table(file=ExClust_Results, sep="\t",header=TRUE,stringsAsFactors=FALSE)
# run the GRangesFromExClustResults function
GRanges.ExClustResults <- GRangesFromExClustResults(results.ExClust=ExClust_Results)
}
\value{
This function returns a GRanges object with ExCluster main function results, in the form of a table of data.
}
\references{
Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118.doi:10.1371/journal.pcbi.1003118
}
\author{
R. Matthew Tanner
}
