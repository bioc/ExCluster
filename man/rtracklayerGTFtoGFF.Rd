\name{rtracklayerGTFtoGFF}
\alias{rtracklayerGTFtoGFF}
\title{
Converts rtracklayer imported GTF files to flattened GFF format
}
\description{
rtracklayerGTFtoGFF converts GTF annotations imported with rtracklayer into flattened GFF files, by converting rtracklayer objects back to GTF files, and then calling the GFF_convert function. ** Please note this function requires 1-2GB of free space to write out a temporary GTF file. This file will be in your tempdir(), and is cleaned up by the function before it finishes. **
}
\usage{
rtracklayerGTFtoGFF(rtracklayer.GTF, GFF.File)
}
\arguments{
    \item{rtracklayer.GTF}{This argument should be given an rtracklayer object which contains a GTF file imported with the tracklayer import() function.
    }
        \item{GFF.File}{A full file path to which the flattened GFF file can be written. This argument is not required, although it is strongly recommended that GFF files be saved for future use. An example file path for a GFF file would be: /Users/username/path/to/file.gff
    }
}
\examples{
# load the sub-sampled GTF file path from the ExCluster package
GTF_file <- system.file("extdata","sub_gen.v23.gtf", package = "ExCluster")
# make an rtrack layer object from this GTF file
rtrack_GTF <- rtracklayer::import(GTF_file)
# now run rtracklayerGTFtoGFF, assigning the results to the GFF object
# in this case we are not specifying GFF.File, and thus not saving the GFF file
GFF <- rtracklayerGTFtoGFF(rtracklayer.GTF=rtrack_GTF)
}
\value{
This function returns a GRanges object with ExCluster main function results, in the form of a table of data.
}
\references{
M. Lawrence, R. Gentleman, V. Carey: "rtracklayer: an {R} package for interfacing with genome browsers". Bioinformatics 25:1841-1842.
}
\author{
R. Matthew Tanner
}
