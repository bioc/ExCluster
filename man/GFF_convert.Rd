\name{GFF_convert}
\alias{GFF_convert}
\title{
Flatten GTF files into GFF format
}
\description{
The GFF_convert function takes a GTF file input and flattens in into a GFF file, which contains discrete, non-overlapping exon bin definitions. Transcripts per exon bin are annotated, and normal gene names are added, in addition the Ensembl IDs.
}
\usage{
GFF_convert(GTF.File,GFF.File)
}
\arguments{
  \item{GTF.File}{A full file path, including file name & extension, to the GTF file which will be flattened. This parameter is required. An example file path for a GTF file would be: /Users/username/path/to/file.gtf
    }
    \item{GFF.File}{A full file path to which the flattened GFF file can be written. This argument is not required, although it is strongly recommended that GFF files be saved for future use. An example file path for a GFF file would be: /Users/username/path/to/file.gff
    }
}
\examples{
# load the sub-sampled GTF file path from the ExCluster package
GTF_file <- system.file("extdata","sub_gen.v23.gtf", package = "ExCluster")
# now run GTF_file without assigning a GFF_file to write out
# assign this output to the GFF variable
GFF <- GFF_convert(GTF.File=GTF_file)
}
\value{
This function returns a GFF data frame, which contains non-overlapping exon bins, flattened from a GTF file. When the function GFF_convert is run, its results should usually be assigned to a variable name, such as 'GFF'.
}
\references{
    Anders, S., Reyes, A., and Huber, W. (2012). Detecting differential usage of exons from RNA-seq data. Genome Research, 22(10):2008-17.
}
\author{
R. Matthew Tanner
}
\note{
The GFF_convert function is based on the ideas from the authors of DEXSeq (Anders et al., 2012). Their dexseq_prepare_annotations.py pipeline for flattening GTF files into GFF format was translated into R. Great thanks are owed to Dylan Siriwardena (MSc) for helping write the initial translated R version.
}
