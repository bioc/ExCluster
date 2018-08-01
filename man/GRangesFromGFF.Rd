\name{GRangesFromGFF}
\alias{GRangesFromGFF}
\title{
Converts GFF formatted data to GRanges format.
}
\description{
GRangesFromGFF converts GFF formatted data generated from the GFF_convert() function into GRanges objects. This data can be formatted either as a dataframe in R, or saved as a GFF file and passed to the function from the GTF.File argument.
}
\usage{
GRangesFromGFF(annot.GFF, GFF.File)
}
\arguments{
    \item{annot.GFF}{If this argument is specified, it must be given a GFF annotation data frame, which is generated from the GFF_convert function. For example, it can be specified as annot.GFF=GFF, assuming your GFF data frame is assigned to the GFF variable. Either annot.GFF or GFF.File must be specified. If both are specified, annot.GFF will take priority.
    }
    \item{GFF.File}{If this argument is specified, it must be given a full file path to a GFF file, including file name and extension. This argument is not required, but either annot.GFF or GFF.File must be specified.
    }
}
\examples{
# grab the path to the sub-sampled example GFF file
GFF_file <- system.file("extdata","sub_gen.v23.ExClust.gff",package="ExCluster")
# run the GRangesFromGFF function
GRanges.GFF <- GRangesFromGFF(GFF.File=GFF_file)
}
\value{
This function returns a GRanges object with GFF annotations, which have been flattened from a GTF file.
}
\references{
Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118.doi:10.1371/journal.pcbi.1003118
}
\author{
R. Matthew Tanner
}
