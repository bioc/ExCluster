ExCluster_errors <- list(

    #################################################################################
    ################################## GFF_convert ##################################
    #################################################################################

    GTF.File_missing =
    "You did not  provide the GFF_convert function a GTF input for an argument.
Please give this function a filepath to a GTF file, such as GTF.File=/path/to/file.gtf,
or provide directly GTF annotations within R to the annot.GTF argument.
These annot.GTF annotations may be a standard data frame or an rtracklayer object",

    bad_GTF_filepath =
        "The GTF file path that you specified does not exist.
Please double check your GTF file's full path & file name.
For example, a GTF.File path for mac/linux users might look like: /User/username/path/to/file.gtf",

    check_GTF_columns =
        "Your GTF file did not have 9 tab delimited columns.
If you have edited your GTF file, please ensure the original number of columns and type of columns are retained.
Additionally, if you have edited your GTF file, please ensure that tabs were used to delimit columns.
If this error persists, please re-download your GTF annotations from the source database, and try again.",

    GTF_missing_features = "Your GTF file did not have 'gene', 'transcript', and 'exon' features present in column 3.
If you have edited your GTF file, please ensure the original number of columns and type of columns are retained.
Also double check that your GTF file was read into R without row or column names, as the GTF format has none.
If this error persists, please re-download your GTF annotations from the source database, and try again.",

    bad_GTF_start_stop = "The start and stop columns (4 & 5) of your GTF file were not properly formatted.
Please ensure your gtf.file was read correctly, without row or column names from the original GTF file.
If this error persists, please re-download your GTF annotations from the source database, and try again.",

    start_stop_not_numeric = "There were negative distances between your stop and start locations.
The start/stop locations of your genomic features in this GTF file appear to be corrupted.
Please ensure column 4 of your GTF file is start location, and column 5 is stop location.
If this error persists, please re-download your GTF annotations from the source database, and try again.",

    bad_GTF_plus_minus = "The strand column of your GTF file (column 7) had more than 2 unique values.
This column should only have '+' and '-' values. Please ensure your GTF file is read into R without row/column names.
If this error persists, please re-download your GTF annotations from the source database, and try again.",

    GTF_plus_minus_missing = "The strand column of your GTF file (column 7) did not contain '+' or '-' values.
This column should only have '+' and/or '-' values. Please ensure your GTF file is read into R without row/column names.
If this error persists, please re-download your GTF annotations from the source database, and try again.",

    GTF_missing_gene_id = "Column 9 in your GTF file did not contain gene_id and/or transcript_id values in the correct places.
Please verify that your GTF file was read into R without row or column names, and that each row has a gene_id in column 9.
If this error persists, please re-download your GTF annotations from the source database, and try again.",

    #################################################################################
    ################################# processCounts #################################
    #################################################################################

    BAM_argument_missing = "No BAM files were entered when running processCounts().
Please ensure the bam.Files argument was assigned a character vector of full BAM file paths.
Alternatively, you can assign the BAM file names as a character vector to the bam.Files variable, such as:
    bam.Files <- c('/path/to/file1.bam','/path/to/file2.bam',...)",

    sample_name_argument_missing = "No sample names were entered to correspond to your BAM files.
Please make sure the sample.Names argument is assigned a character vector when calling processCounts().
Alternatively, you can assign the character vector to the sample.Names variable as follows:
sample.Names <- c('control1_rep1','control1_rep2',...)",

    paired_reads_not_logical = "You did not correctly enter the 'pairedReads' variable as TRUE or FALSE.
Please use pairedReads=TRUE or pairedReads=FALSE when calling processCounts()",

    BAM_sample_length_unequal = "The number of BAM files and sample names are not the same!
Please ensure the lengths of bam.Files and sample.Names variables are the same.",

    GFF_file_doesnt_exist = "The GFF file path you provided to GFF.File does not exist.
Please verify you have given the correct, full filepath including file extension.
    For example, a file path may look like: /Users/username/path/to/file.gff",

    GFF_file_argument_missing = "The GFF.File argument was not provided when processCounts was run, nor was the annot.GFF argument assigned.
One of annot.GFF or GFF.File must be specified when running processCounts().
The GFF.File argument should be a full file path to a GFF file, such as: /Users/username/path/to/file.gff
Alternatively, you can give the annot.GFF argument a GFF data frame from the GFF_convert function, such as: annot.GFF = GFF",

    GFF_missing_columns = "Your GFF file did not have 9 columns.
Please re-run GFF_convert from an unaltered GTF file, and pass that output to processCounts.
Consult the vignette if you are unsure how to proceed.",

    GFF_missing_GFF3_fields = "Your GFF file did not have ID, Name, and Transcript fields properly formatted in column 9.
Please ensure you are using a GFF file produced from the ExCluster GFF_convert function.
Unfortunately this function can only accept GFF3 formatted GFF files from GFF_convert.",

    non_integer_numcores = "You entered a non-integer value for num.Cores -- it should have a value of >= 1.
Defaulting number of cores to 1 for this run.",

    temp_dir_inaccessible = "The temp.Dir path you provided was not accessible by R for read/write.
    Please ensure the temp.Dir argument is given a valid folder path with read/write permissions.
    Defaulting temporary directory to R's tempdir() function.",

    low_exon_bin_count = "You had fewer than 1000 exon bins with RNA-seq reads counted into them across conditions.
Please be aware that ExCluster expects whole transcriptome GFF and BAM files counted, not individual chromosome data.
Without these full files, ExCluster cannot properly correct library sizes and tune statistics.
    >> If you are simply testing an ExCluster example, please ignore this warning.",

    #################################################################################
    ################################### ExCluster ###################################
    #################################################################################

    exon_counts_missing = "You did provide a normalized exon count matrix to the exonCounts argument.
Please obtain normalized exon counts from processCounts() and assign the resutls to exonCounts.
    For example, if you followed the vignette, you should have a normCounts variable.
    In this case, you could specify: exonCounts=normCounts",

    cond_nums_missing = "You did not provide condition numbers to the cond.Nums argument.
Please ensure you specify which columns of your exon read counts correspond to which condition.
For example, specify cond.Nums=c(1,1,1,2,2,2) if you have 6 samples between 2 conditions.
Please note the order of the numbers must correspond to the order of the columns.",

    exon_count_length_incorrect =
"Your number of exonCounts columns did not equal the length of your number of condition identifiers.
Please ensure that you entered correct condition identifiers to cond.Nums, such as cond.Nums=c(1,1,1,2,2,2).
Additionally, please also ensure your normalized counts provided to exonCounts have only count data columns.
All Gene:ExonBin identifiers should be provided as rownames. If you are unsure, please refer to the vignette.",

    exon_read_counts_duplicated = "One of your normalized exon count columns was duplicated!
If this is a mistake, please double check your normalized count table, or re-run processCounts().
However, if you were attempting to duplicate columns to analyze data without biological replicates, this will not work.
ExClust requires biological replicates to estimate variance within conditions, and duplicate columns yield variances of zero.",

    GFF_annotations_missing = "You must provide either GFF annotation data or a GFF file path.
Please assign the annot.GFF argument a variable name containing the GFF data, such as annot.GFF=GFF.
Alternatively, please specify a full GFF file path, such as: GFF.File=/Users/username/path/to/file.gff",

    GFF_file_inaccessible = "The GFF file path you entered did not exist.
Please verify that you have assigned the exactly correct path to the GFF.File argument.
For example, your argument when running ExCluster should look something like this: GFF.File=/Users/username/path/to/file.gff",

    improper_cond_nums = "You did not specify exactly 2 condition numbers to the cond.Nums argument.
For example, if your data contains two conditions with three replicates each, please use: cond.Nums=c(1,1,1,2,2,2).
If your count data alternates condition columns, such as cond1, cond2, cond1, cond2, etc., use: cond.Nums=c(1,2,1,2,1,2)",

    not_enough_replicates_cond1 = "You only provided biological replicate for your first condition.
ExCluster, unfortunately, requires at least 2 biological replicates per condition to function.
If this was a mistake, please re-enter your cond.Nums argument with at least 2 replicates per condition number.",

    not_enough_replicates_cond2 = "You only provided biological replicate for your second condition.
ExCluster, unfortunately, requires at least 2 biological replicates per condition to function.
If this was a mistake, please re-enter your cond.Nums argument with at least 2 replicates per condition number.",

    ExClust_GFF_bad_colnum = "You did not enter a valid GFF annotation data frame or a proper path to your GFF file.
It is possible that your GFF file is not properly formatted with the expected 9 columns.
Please re-run GFF_convert and then re-assign either annot.GFF or GFF.File when calling ExCluster.
See the ExCluster manual and vignette for more information.",

    GFF_missing_exon_feature = "Your GFF file did not have column 3 beginning with 'exon' in all cases.
Please ensure your GFF file is produced from the GFF_convert function, and read into R without row/column names.
If this error persists, please re-produce your GFF file from a GTF file using GFF_convert.",

    GFF_numeric_error = "The GFF file start and end genomic coordinate columns (4 & 5) were not properly formatted numbers.
Please ensure your GFF annotations were produced from the GFF_convert function.
If you are reading saved GFF annotations from a file, please ensure no row or column names were applied to the file.
If this error persists, please re-produce your GFF file from a GTF file using GFF_convert.",

    GFF_negative_dists = "There were negative distances between your stop and start locations in your GFF file.
Please ensure your GFF annotations were produced from the GFF_convert function.
If you are reading saved GFF annotations from a file, please ensure no row or column names were applied to the file.
Columns 4 & 5 of the GFF file format should contain start/stop locations with only numeric values throughout.
If this error persists, please re-produce your GFF file from a GTF file using GFF_convert.",

    bad_GFF_strand_column = "The strand column of your GFF file (column 7) had more than 2 unique values.
This column should only have '+' and '-' values. Please ensure your GFF file is read into R without row/column names.
If this error persists, please re-produce your GFF file from a GTF file using GFF_convert.",

    GFF_plus_minus_missing = "The strand column of your GTF file (column 7) did not contain '+' or '-' values.
This column should only have '+' and '-' values. Please ensure your GFF file is read into R without row/column names.
If this error persists, please re-produce your GFF file from a GTF file using GFF_convert.",

    match_GFF_counts_failed = "The GFF file you provided does not match the GFF file used to count exonic reads.
Please consider converting your GTF file to GFF again, and then re-counting exon reads.
Otherwise, please make sure the same GFF file used to count exon reads is provided to ExCluster's input arguments.",

    too_few_pvalues = "You have entered transcriptome data which contains fewer than 5000 expressed genes.
ExCluster cannot run advanced statistics (yielding better tuned FDRs) with < 5000 expressed genes.
Please avoid running ExCluster on sub-sections of a genome (i.e. only certain chromosomes),
as background expressed genes are used to tune statistics.

If you are running test data to get ExCluster working, ignore this warning.",

    ExClust_bad_outDir = "You specified an out.Dir directory for ExCluster that could not be written to.
    Please provide a valid directory that ExCluster can write files to.
    ExCluster will finish running, however your file has not been written.",

    FDR_cutoff_too_high = "FDR.cutoff was assigned higher than 0.2, which is not allowed.
    Using FDR cutoffs above 20% generates unnecessarily high false discovery rates.
    Defaulting to 0.05 FDR for this run of ExCluster().",

    plot_results_no_outDir = paste("You specified that plot.Results=TRUE, ",
                "but out.Dir was not specified or could not be written to. ",
                "out.Dir must be assigned a valid directory that can be created and/or written to. ",
                "Your exon log2FCs have not been plotted, however, you may re-plot them ",
                "with the plotExClustResults function -- this requires saved ExClust results. ",
                "Please consult the vignette for more instructions.", sep=""),

    plot_type_failure = paste("ExCluster was unable to plot either PNG or bitmap. ",
                              "This error commonly occurs because your operating system ",
                              "lacks both Ghostscript and X11 forwarding. Please ensure ",
                              "R has access to either Ghostscrpit or X11 forwarding.",sep="")

)
