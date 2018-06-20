\name{processCounts}
\alias{processCounts}
\title{
Counts reads per exon bin for each BAM file, using GFF exon annotations -- also normalizes library sizes.
}
\description{
The processCounts function takes an input of BAM files, and counts reads per exon bin from GFF annotations. GFF annotations may be given either in the form of a data frame assigned to annot.GFF, or the path to a GFF file given to GFF.File. processCounts also requires sample names per BAM file, to unambiguously name the read count columns. Reads may be counted as pairs or individually for paired-end reads, by setting the pairedReads parameter to TRUE or FALSE.
}
\usage{
processCounts(bam.Files,sample.Names,annot.GFF,GFF.File,pairedReads,out.File,temp.Dir=NULL,num.Cores=NULL)
}
\arguments{
  \item{bam.Files}{An array of full file paths, specifying BAM files from RNA-seq experiments for which exon bin reads will be counted. Each item in the array must specify a full file path to a unique BAM file -- duplicate BAM files are not allowed. This is a required input.
    }
    \item{sample.Names}{An array of sample names corresponding, in exact order, to each BAM file provided to the bam.Files argument. The length of sample.Names must also match the length of bam.Files exactly. This is a required input.
    }
    \item{annot.GFF}{If this argument is specified, it must be given a GFF annotation data frame, which is generated from the GFF_convert function. For example, it can be specified as annot.GFF=GFF, assuming your GFF data frame is assigned to the GFF variable. Either annot.GFF or GFF.File must be specified. If both are specified, annot.GFF will take priority.
    }
    \item{GFF.File}{If this argument is specified, it must be given a full file path to a GFF file, including file name and extension. For example, this may look like: /Users/username/path/to/file.gff This argument is not required, but either annot.GFF or GFF.File must be specified.
    }
    \item{pairedReads}{pairedReads should be assigned a value of either TRUE or FALSE, and functions as a logical operator to determine whether or not reads should be counted as pairs. For example, if you are counting reads from a BAM file with paired-end reads, you should set pairedReads=TRUE By default pairedReads is FALSE, and thus need not be set for single-end reads.
    }
    \item{out.File}{This argument should be given a variable with a string which contains the full file path to which you wish to write your normalized exon bin counts. For example, you could set outFile<-"/path/to/normCounts.txt", and run processCounts with out.File=outFile
    }
    \item{temp.Dir}{This argument should be given a variable with a string which contains the full file path to a temporary directory to which you wish to write featureCounts temporary files. In some cases this can greatly speed up read counting. By default, featureCounts temporary files will be written to R's tempdir() directory.
    }
    \item{num.Cores}{This argument should be given a positive integer value which corresponds to the number of cores you with featureCounts to use, when exon reads are counted from your GFF file. For example, using 4 cores would look like: num.Cores=4  This argument defaults to using 1 core within the function.}
    }
}
\examples{
# specify the path to the ExCluster package
ExClust_Path <- system.file(package="ExCluster")
# now find the bam files within that folder
bamFiles <- list.files(ExClust_Path,recursive=TRUE,pattern="*.bam",
    full.names=TRUE)
# now grab the path to the sub-sampled example GFF file
GFF_file <- system.file("extdata","sub_gen.v23.ExClust.gff",
    package="ExCluster")
# assign sample names (only 2 replicates per condition in this example)
sampleNames <- c("iPSC_cond1_rep1","iPSC_cond1_rep2","iPSC_cond2_rep1",
    "iPSC_cond2_rep2")
# now run processCounts, with pairedReads=TRUE for paired-end data
normCounts <- processCounts(bam.Files=bamFiles, sample.Names=sampleNames,
    GFF.File=GFF_file, pairedReads=TRUE)
}
\value{
This function returns normalized exon counts per non-overlapping exon bin in the form of a data frame, in which the user has manually entered sample names which will be present as column headers. This should typically be assigned to a variable when run, such as 'normCounts'.
}
\references{
Liao, Y., Smyth, G.K., and Shi, W. (2013). The Subread aligner: fast, accurate and scalable read
mapping by seed-and-vote. Nucleic Acids Research, 41(10):e108.
}
\author{
R. Matthew Tanner
}
\note{
processCounts runs the featureCounts function from the Rsubread package (Liao et al., 2013). As their function is being called from an external package, the Rsubread authors and license owners make no warranties as to its performance.
}
