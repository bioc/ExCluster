plotExonlog2FC <- function(results.Data=NULL, out.Dir=NULL, FDR.cutoff=0.05, plot.Type="bitmap"){

    ### make sure R doesn't add factors -- R creators never should have made this default = TRUE
    options(stringsAsFactors=FALSE)
    # check to make sure results.Data parameter is given as input
    if (is.null(results.Data) == TRUE){
        stop(call="You must assign the results.Data argument a variable containing the results of an ExCluster analysis.
    For example, if previous ExCluster results were assigned to clustResults, you would specify: results.Data=clustResults")
    }
    # check to make sure out.Dir is specified
    if (is.null(out.Dir) == TRUE){
        stop(call="plotExClustResults requires a direcotry to write images to, but was not specified.
    Please specify the full file path to a folder to which exon log2FC plots may be written.")
    }else{
        # now that out.Dir is specified, make sure its parent directory exists
        if (dir.exists(out.Dir) == FALSE){
            # try to make out.Dir no matter what
            try(dir.create(out.Dir, recursive=TRUE, showWarnings = FALSE), silent=TRUE)
        }
    }

    ### now make sure we can write to the out.Dir
    WriteCheck <- file.access(out.Dir, mode=2)
    if (WriteCheck != 0){
        stop(call="ExCluster did not have write priveledges for the out.Dir you specified.
             Please ensure you specify a valid directory that can be created and written to by R.")
    }

    ### now check bitmap and PNG settings for plot.Type
    if (is.null(plot.Type) == TRUE){
        plot.Type <- "bitmap"
    }else if (plot.Type != "bitmap" && plot.Type != "PNG"){
        plot.Type <- "bitmap"
    }

    ### now check that PNG and/or bitmap can be plotted
    # plot.Check is FALSE by default, unless proven TRUE
    plot.Check <- FALSE
    if (plot.Type == "bitmap"){
        # check if bitmap can be plotted
        test.bitmap <- substr(testBMplot(out.Dir)[1],1,5)
        if (test.bitmap == "Error"){
            # check if PNG can be plotted
            test.PNG <- substr(testPNGplot(out.Dir)[1],1,5)
            plot.Type <- "PNG"
            if (test.PNG != "Error"){
                plot.Check <- TRUE
            }
        }else{
            plot.Check <- TRUE
        }
    }else{
        # if plot.Type == "PNG", repeat the steps above in reverse
        test.PNG <- substr(testPNGplot(out.Dir)[1],1,5)
        if (test.PNG == "Error"){
            test.bitmap <- substr(testBMplot(out.Dir)[1],1,5)
            plot.Type <- "bitmap"
            if (test.bitmap != "Error"){
                plot.Check <- TRUE
            }
        }else{
            plot.Check <- TRUE
        }
    }

    ### make sure either bitmap or PNG passed
    if (plot.Check == FALSE){
        stop(call = ExCluster_errors$plot_type_failure)
    }

    ############################# FUNCTIONS ###############################

    ErrorBar <- function(x,y1,y2,MaxXvalue,BinSize){
        ErrorWidth <- min((MaxXvalue*0.015),(BinSize/3))
        lines(c(x,x),c(y1,y2))
        lines(c(x-ErrorWidth,x+ErrorWidth),c(y1,y1))
        lines(c(x-ErrorWidth,x+ErrorWidth),c(y2,y2))
    }

    ########################### NOW PREPARE DATA ###########################

    ### sort results.Data
    results.Data <- results.Data[order(results.Data$EnsG),]

    ### Target genes of interest
    TargetGenes <- unique(results.Data$EnsG[which(results.Data$FDR <= FDR.cutoff)])

    if (length(TargetGenes) > 0){
        ### grab IDs to make data parsing faster
        IDs <- gsub("\\:.*","",results.Data$EnsG)
        IDs <- gsub("\\..*","",IDs)
        IDs <- substr(IDs,(nchar(IDs[1])-9),nchar(IDs[1]))
        # Now add the IDs to results.Data
        results.Data$IDs <- IDs

        ### Do the same for TargetGenes
        IDs <- gsub("\\:.*","",TargetGenes)
        IDs <- gsub("\\..*","",IDs)
        IDs <- substr(IDs,(nchar(IDs[1])-9),nchar(IDs[1]))

        ### loop through results looking for gene data so we can grab FDR, min/max results.Data, etc.
        Start <- 1
        Counter <- 1
        RemovalIndices <- NULL

        Start <- 1

        for (i in 2:nrow(results.Data)){
            # check to see if we have reached a new gene
            if ((results.Data$IDs[i] != results.Data$IDs[(i-1)]) == TRUE | i == nrow(results.Data)){
                # determine stop of gene
                if (i == nrow(results.Data)){
                    Stop <- i
                }else{
                    Stop <- (i-1)
                }
                # now check to see if the current gene of interest is in TargetGenes
                if (results.Data$IDs[Start] == IDs[Counter]){
                    Counter <- min((Counter + 1),length(IDs))
                }else{
                    RemovalIndices <- c(RemovalIndices, c(Start:Stop))
                }
                # reset start for next gene
                Start <- i
            }
        }

        if (length(RemovalIndices) > 0){
            ### Now remove the rows specified in RemovalIndices to yield significant gene data
            sig.Results <- results.Data[-c(RemovalIndices),]
            rm(results.Data)
        }else{
            sig.Results <- results.Data
        }
        # calculate length of each exon bin
        sig.Results$End <- as.numeric(sig.Results$End)
        sig.Results$Start <- as.numeric(sig.Results$Start)
        sig.Results$Length <- sig.Results$End - sig.Results$Start + 1

        ############################# NOW PLOT DATA #############################

        ### loop through signficant results, plotting each gene
        Start <- 1
        for (i in 2:nrow(sig.Results)){
            # check to see if we have reached a new gene
            if ((sig.Results$IDs[i] != sig.Results$IDs[(i-1)]) == TRUE | i == nrow(sig.Results)){
                # determine stop of gene
                if (i == nrow(sig.Results)){
                    Stop <- i
                }else{
                    Stop <- (i-1)
                }
                ### Grab temporary gene data to plot
                TempData <- sig.Results[Start:Stop,]
                # Set up data variables to be used in plotting
                TargetGene <- TempData$EnsG[1]
                TargetName <- TempData$Gene_name[1]
                TargetFDR <- TempData$FDR[1]
                # ceiling of abs(log2FC) + SD + 1 log2FC to set for +/- y-axis limits
                Maxlog2FC <- ceiling(max(abs(TempData$log2FC) + sqrt(TempData$log2Var)))
                # maximum x-axis
                MaxX <- sum(TempData$Length) + nrow(TempData) - 1

                ### now set up png for plotting
                # check to see if we are plotting bitmap or PNG
                if (plot.Type == "bitmap"){
                    bitmap(file=paste(out.Dir,"/",TargetGene,"_",TargetName,".png",sep=""), type="png16m",
                           width=6,height=4,units="in",res = 400, taa=4,gaa=4)
                }
                if (plot.Type == "PNG"){
                    png(filename = paste(out.Dir,"/",TargetGene,"_",TargetName,".png",sep=""),
                        width=6,height=4,units="in",res = 400)
                }
                par(ps = 12, cex = 1, cex.main = 1,mar=(c(5, 5, 4, 2)+0.2))
                plot(c(0,MaxX), c(-(Maxlog2FC),Maxlog2FC) , type= "n",
                     ylab = "log2 fold change (log2FC)",xlab = "nucleotide position (number of base pairs)",
                     main = paste(TargetName,"   ",TargetGene,"    ","FDR=",TargetFDR,sep=""),cex.main=0.9)
                # add a light grey line at log2FC = 0
                lines(x = c(-500,MaxX+500), y = c(0,0), col = "grey")
                # run plotting code
                # Start of Exon = 0
                ExonStart <- 0
                # Bar height
                BH <- Maxlog2FC/60
                for (j in seq(nrow(TempData))){
                    rect(ExonStart, TempData$log2FC[j]-BH, ExonStart + TempData$Length[j],
                         TempData$log2FC[j]+BH, col = "red")
                    ErrorBar((ExonStart + (TempData$Length[j] / 2)),(TempData$log2FC[j]-sqrt(TempData$log2Var[j])),
                             (TempData$log2FC[j]+sqrt(TempData$log2Var[j])),MaxX,TempData$Length[j])
                    ExonStart <- ExonStart + TempData$Length[j] + 1
                }
                dev.off()
                # reset the Start variable for the next gene
                Start <- i
            }
        }
    }
}

