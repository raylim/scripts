#!/usr/bin/env Rscript
# runs HMM copy on tumor normal wig files

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("HMMcopy"));
suppressPackageStartupMessages(library("hwriter"));

optionList <- list(make_option(c("-g", "--gcWig"), default = NULL, help = "GC content wig file"),
                   make_option(c("-m", "--mapWig"), default = NULL, help = "Mappability wig file"),
                   make_option(c("-n", "--normalWig"), default = NULL, help = "Normal wig file"),
                   make_option(c("-o", "--out"), default = NULL, help = "Output file"),
                   make_option(c("-d", "--outDir"), default = NULL, help = "Output directory"));
posArgs <- c('tumor_wig')
parser <- OptionParser(usage = paste('%prog [options]', paste(posArgs, collapse=' ')),  option_list=optionList)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args) != 1) {
	print_help(parser)
	print(arguments$args)
	stop('Incorrect number of required positional arguments')
} else if (is.null(opt$gcWig)) {
    print_help(parser)
    print(arguments$args)
    stop('Need GC content wig file')
} else if (is.null(opt$out)) {
    print_help(parser)
    print(arguments$args)
    stop('Need output file')
} else if (is.null(opt$mapWig)) {
    print_help(parser)
    print(arguments$args)
    stop('Need mappability wig file')
}

cmdArgs <- arguments$args
for (i in 1:length(cmdArgs)){
    assign(posArgs[i], cmdArgs[i])
}


tumorWigFile <- cmdArgs[1]
normalWigFile <- opt$normalWig
gcWigFile <- opt$gcWig
mapWigFile <- opt$mapWig

tumorData <- wigsToRangedData(tumorWigFile, gcWigFile, mapWigFile)
tumorData <- correctReadcount(tumorData)

if (!is.null(normalWigFile)) {
    normalData <- wigsToRangedData(normalWigFile, gcWigFile, mapWigFile)
    normalData <- correctReadcount(normalData)
    tumorData$copy <- tumorData$copy - normalData$copy
}

param <- HMMsegment(tumorData, getparam = TRUE)
param$mu <- log(c(1, 1.4, 2, 2.7, 3, 4.5) / 2, 2) # Gavin's parameters
param$m <- param$mu

tumorSeg <- HMMsegment(tumorData, param)
write.table(as.data.frame(tumorSeg$segs), file = opt$out, sep = '\t', quote = F)

if (!is.null(opt$outDir)) {
    dir.create(opt$outDir, showWarnings = F)
    setwd(opt$outDir)
    dir.create('graphics', showWarnings = F)
    pg <- openPage('index.html')

    gfn <- 'graphics/seg_plot.png'
    png(gfn, height = 800, width = 1000)
    plotSegments(tumorData, tumorSeg, pch = '.', ylab = "Tumour Copy Number", xlab = "Chromosome Position")
    null <- dev.off()
    hwriteImage(gfn, pg)

    gfn <- 'graphics/bias_plot.png'
    png(gfn, height = 1000, width = 1000)
    plotBias(tumorData, pch = 19)
    null <- dev.off()
    hwriteImage(gfn, pg)

    for (chr in paste("chr", 1:21, sep = '')) {
        gfn <- paste('graphics/seg_plot_', chr, '.png', sep = '')
        png(gfn, height = 1000, width = 1000)
        plotSegments(tumorData, tumorSeg, pch = '.', chr = chr, ylab = "Tumour Copy Number", xlab = "Chromosome Position", main = chr)
        cols <- stateCols()
        legend("topleft", c("HOMD", "HETD", "NEUT", "GAIN", "AMPL", "HLAMP"), fill = cols, horiz = T, bty = "n", cex = 0.5)
        null <- dev.off()
        hwriteImage(gfn, pg)
    }

    for (chr in paste("chr", 1:21, sep = '')) {
        gfn <- paste('graphics/cor_plot_', chr, '.png', sep = '')
        png(gfn, height = 1000, width = 1000)
        plotCorrection(tumorData, pch = '.', chr = chr)
        null <- dev.off()
        hwriteImage(gfn, pg)
    }

    closePage(pg)
}

