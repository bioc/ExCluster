ExClust_compute_stats <- function(gene.IDs=NULL, gene.PVals=NULL){

    ### first we need to grab the gene names & p-values for each gene
    # grab the 1st row of each gene
    pvalue.Indices <- match(unique(gene.IDs), gene.IDs)
    # now grab the gene IDs from these rows
    NamesGenePVals <- gene.IDs[pvalue.Indices]
    # lastly grab the p-values from these rows
    GenePVals <- gene.PVals[pvalue.Indices]

    ### ExCluster progression message
    message("Running statistical tests ...","",sep="\n")

    ### data frame for testing p-values
    Catch_PVals <- data.frame(NamesGenePVals,GenePVals)
    colnames(Catch_PVals) <- c("EnsG","pvals")
    # remove NAs
    Catch_PVals <- Catch_PVals[which(complete.cases(Catch_PVals[,2])),]
    # sort PVals from lowest to highest
    Catch_PVals <- Catch_PVals[order(Catch_PVals[,2]),]

    if (nrow(Catch_PVals) >= 5000){
        ### normalize p-values
        # estimate the number of 'truly spliced genes'
        pvalCutoff <- EstNullHypo(Catch_PVals[,2])
        # estimate the number of truly spliced genes (may include some false positives)
        NumTrue <- length(which(Catch_PVals[,2] <= pvalCutoff))
        # estimate the number of null hypothesis genes
        NumNull <- length(which(Catch_PVals[,2] > pvalCutoff))
        # the estimated number of false positives in 'NumTrue'
        NumFalsePos <- (length(which(Catch_PVals[,2] > pvalCutoff & Catch_PVals[,2] <= 0.1)))/
                        ((0.10-pvalCutoff)/pvalCutoff)
        # now estimate the minimum p-value of the null hypothesis
        minNullPVal <- gm_mean(vapply(seq(100),function(x){
            # generate an initial uniform p-value distribution for the total NumNull
            InitialNullPVals <- sort(runif((NumNull+NumFalsePos),0,1))
            # now use this to estimate the minimum p-value of the null hypothesis distribution
            return(InitialNullPVals[NumFalsePos+1])
        },
        FUN.VALUE = c(Res=0)))

        ### if we have > 0 NumTrue p-values
        if (NumTrue > 0){
            # Now give the null hypothesis p-values uniform distriubtions
            Catch_PVals[seq((NumTrue+1),(NumTrue+NumNull)),2] <- sort(runif(n=NumNull,0,1))
            # Now divide all true p-values by a factor based on the lowest null hypothesis value
            Catch_PVals[seq(NumTrue),2] <- Catch_PVals[seq(NumTrue),2]/
                                            (Catch_PVals[NumTrue,2]/Catch_PVals[(NumTrue+1),2])
            # Now compute FDR values
            FDRs <- FDRcalc(Catch_PVals[,2],NumTrue,NumNull)
        }else{
            ### if we have no truly significant p-values from our estimate
            Catch_PVals[,2] <- sort(runif(n=NumNull,0,1))
            FDRs <- rep(1,nrow(Catch_PVals))
        }

    }else{

        ### This part of the code will be run if less than 5000 expressed genes were identified
        warning(call = ExCluster_errors$too_few_pvalues)

        # because we have too few pvalues (< 5000), we cannot run our advanced stats
        # we therefore run a basic correction with Benjamini-Hochberg
        pvalCutoff <- 0.025
        # estimate the number of truly spliced genes
        NumTrue <- nrow(Catch_PVals[Catch_PVals[,2] <= pvalCutoff,])
        # estimate the number of null hypothesis genes
        NumNull <- nrow(Catch_PVals[Catch_PVals[,2] > pvalCutoff,])
        # set the minimum value a null hypothesis p-value can be corrected to
        minNullPVal <- 0.005

        ### if we have truly differentially spliced genes
        if (NumTrue > 0){
            # Now give the null hypothesis p-values uniform distriubtions
            Catch_PVals[seq((NumTrue+1),(NumTrue+NumNull)),2] <- sort(runif(n=NumNull,minNullPVal,1))
            # Now divide all true p-values by a factor based on the lowest null hypothesis value
            Catch_PVals[seq(NumTrue),2] <- Catch_PVals[seq(NumTrue),2]/(Catch_PVals[NumTrue,2]/minNullPVal)
            # Now compute FDR values
            FDRs <- p.adjust(p = Catch_PVals[,2],method="BH")
        }else{
            ### if we have no truly significant p-values from our estimate
            Catch_PVals[,2] <- sort(runif(n=NumNull,0,1))
            FDRs <- rep(1,nrow(Catch_PVals))
        }
    }

    ### Make statistics table per gene with pvals and FDRs
    Stats_table <- data.frame(Catch_PVals,FDRs)
    return(Stats_table)
}
