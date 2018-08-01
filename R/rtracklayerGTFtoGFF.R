rtracklayerGTFtoGFF <- function(rtracklayer.GTF=NULL, GFF.File=NULL){

    ### This function converts rtracklayer imported GTF files into flattened GFF files
    ### it functions by converting rtracklayer format back to GTF, and then calls GFF_convert()

    ### export rtracklayer object to GTF file (temporary gzip file)
    export(rtracklayer.GTF, format = "GTF", con=paste(tempdir(),"/temp.gtf",sep=""))

    ### now read this file back in
    annot.GTF <- read.table(file=paste(tempdir(),"/temp.gtf",sep=""),sep="\t", stringsAsFactors=FALSE)

    ### now collapse exons into GFF annotations
    annot.GFF <- GFF_convert(GTF.File=paste(tempdir(),"/temp.gtf",sep=""), GFF.File=GFF.File)

    ### clean up
    file.remove(paste(tempdir(),"/temp.gtf",sep=""))

    ### return result & end function
    return(annot.GFF)
}
