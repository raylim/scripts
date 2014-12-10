#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("TitanCNA"));
suppressPackageStartupMessages(library("rtracklayer"));
suppressPackageStartupMessages(library("doMC"));
suppressPackageStartupMessages(library("hwriter"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
        make_option("--outPrefix", default = NULL, type = "character", action = "store", help ="output prefix (required)"),
        make_option("--plotPrefix", default = NULL, type = "character", action = "store", help ="plot output prefix (required)"),
        make_option("--gcWig", default = NULL, type = "character", action = "store", help ="GC wig (required)"),
        make_option("--mapWig", default = NULL, type = "character", action = "store", help ="mappability wig (required)"),
        make_option("--numCores", default = 1, type = "integer", action = "store", help ="number of cores [default = %default]"),
        make_option("--numClusters", default = 5, type = "integer", action = "store", help ="number of clusters [default = %default]"),
        make_option("--tumorWig", default = NULL, type = "character", action = "store", help ="tumor wig (required)"),
        make_option("--normalWig", default = NULL, type = "character", action = "store", help ="normal wig (required)"),
        make_option("--includeY", default = F, action = "store_true", help ="include Y chromosome"),
        make_option("--targetBed", default = NULL, type = "character", action = "store", help ="targeted interval bed"))

parser <- OptionParser(usage = "%prog [options] [tumour allele count file]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need tumour allele count file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$plotPrefix)) {
    cat("Need plot output prefix\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$gcWig)) {
    cat("Need gc wig file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$mapWig)) {
    cat("Need mappability wig file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$tumorWig)) {
    cat("Need tumour wig file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$normalWig)) {
    cat("Need normal wig file\n\n")
    print_help(parser);
    stop();
}

#options(cores = opt$numCores)
registerDoMC(opt$numCores)

chroms <- c(1:22, "X")
if (opt$includeY) {
    chroms <- c(chroms, "Y")
}

pg <- openPage(paste(opt$outPrefix, '_titan_report_', opt$numClusters, '.html', sep = ''), title = 'TITAN Plots')

fn <- arguments$args[1]
Data <- loadAlleleCounts(fn, header = F)
params <- loadDefaultParameters(copyNumber=5, numberClonalClusters=opt$numClusters, symmetric=TRUE, data = Data)


if (!is.null(opt$targetBed)) {
    targetGr <- import(opt$targetBed)
    targets <- as.data.frame(targetGr)[,1:3]
    colnames(targets) <- c("chr", "start", "stop")
    cnData <- correctReadDepth(opt$tumorWig, opt$normalWig, opt$gcWig, opt$mapWig, targetedSequence = targets)
} else {
    cnData <- correctReadDepth(opt$tumorWig, opt$normalWig, opt$gcWig, opt$mapWig)
}

logR <- getPositionOverlap(Data$chr, Data$posn, cnData)
Data$logR <- log(2^logR)
rm(logR, cnData)

Data <- filterData(Data, chroms, minDepth = 10, maxDepth = 250)
mScore <- as.data.frame(wigToRangedData(opt$mapWig))
mScore <- getPositionOverlap(Data$chr, Data$posn, mScore[,-4])
Data <- filterData(Data, chroms, minDepth = 10, maxDepth = 250, map = mScore, mapThres = 0.8)

convergeParams <- runEMclonalCN(Data, gParams=params$genotypeParams, nParams=params$normalParams,
                                pParams=params$ploidyParams, sParams=params$cellPrevParams,
                                maxiter=20, maxiterUpdate=1500, txnExpLen=1e15, txnZstrength=1e5,
                                useOutlierState=FALSE,
                                normalEstimateMethod="map", estimateS=TRUE, estimatePloidy=TRUE)


optimalPath <- viterbiClonalCN(Data, convergeParams)

fn <- paste(opt$outPrefix, '.titan_', opt$numClusters, ".txt", sep = "")
if (opt$numClusters <= 2) {
    results <- outputTitanResults(Data, convergeParams, optimalPath, filename = fn, posteriorProbs = F, subcloneProfiles = T)
} else {
    results <- outputTitanResults(Data, convergeParams, optimalPath, filename = fn, posteriorProbs = F)
}

fn <- paste(opt$outPrefix, '.params_', opt$numClusters, ".txt", sep = "")
outputModelParameters(convergeParams, results, fn)

# plots
norm <- convergeParams$n[length(convergeParams$n)]
ploidy <- convergeParams$phi[length(convergeParams$phi)]

#library(SNPchip)  ## use this library to plot chromosome idiogram (optional)
for (chr in intersect(results$Chr, chroms)) {
    outplot <- paste(opt$plotPrefix, '.titan_', opt$numClusters, ".chr", chr, ".png", sep = '')
    hwriteImage(basename(output), pg, br = T)
    png(outplot,width=1200,height=1000,res=100, type = 'cairo-png')
    if (opt$numClusters <= 2) { 
        par(mfrow=c(4,1))
    } else {
        par(mfrow=c(3,1))
    }
    plotCNlogRByChr(results, chr, ploidy=ploidy, geneAnnot=NULL, spacing=4,ylim=c(-4,6),cex=0.5,main= paste("Chr", chr))
    plotAllelicRatio(results, chr, geneAnnot=NULL, spacing=4, ylim=c(0,1),cex=0.5,main=paste("chr", chr))
    plotClonalFrequency(results, chr, normal=tail(convergeParams$n,1), geneAnnot=NULL, spacing=4,ylim=c(0,1),cex=0.5,main= paste("Chr", chr))
    if (opt$numClusters <= 2){ 
        plotSubcloneProfiles(results, chr, cex = 2, spacing=6, main=paste("Chr", chr))
    }
    #pI <- plotIdiogram(chr,build="hg19",unit="bp",label.y=-4.25,new=FALSE,ylim=c(-2,-1))
    null <- dev.off()
}

closePage(pg)

