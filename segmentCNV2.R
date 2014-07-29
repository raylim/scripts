#!/usr/bin/env Rscript
# segment copy number data and generate plot

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("CGHcall"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--centromereFile", default = NULL, type = "character", action = "store", help ="centromere file"),
                make_option("--prefix", default = NULL, type = "character", action = "store", help ="Output prefix (required)"))

parser <- OptionParser(usage = "%prog [options] inDir", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 1) {
    cat("Need copy number file\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$prefix)) {
    cat("Need output prefix\n\n");
    print_help(parser);
    stop();
} else {
    cnFile <- arguments$args[1];
}

cn <- read.table(cnFile, header=T, as.is=T)
keep <- which(cn[,1] %in% c(1:22, "X"))
if (length(rm) > 0) { cn <- cn[keep,]}
cn[which(cn[,1]=="X"),1] <- 23
cn[,1] <- as.numeric(cn[,1])
cn <- cn[order(cn[,1], cn[,2]),]
cn <- cbind(name = paste(cn[,1], cn[,2], sep="_"), cn[,c(1:3,7)])
cgh <- make_cghRaw(cn)
normalized <- normalize(cgh, smoothOutliers=T, trim=0.025, smooth.region=10, outlier.SD.scale=2.5)
segmented <- segmentData(normalized, relSDlong=3, undo.splits="sdundo", undo.SD=2, alpha=0.000001, trim=0.025)

fn <- paste(opt$prefix, '.segment.Rdata', sep = '')
save(segmented, file = fn)
Data <- cbind(fData(segmented), copynumber(segmented), segmented(segmented))
colnames(Data)[5] <- "log2_ratio_seg"
write.table(Data, file = paste(opt$prefix, ".seg.txt", sep=""), col.names=NA, quote=F, sep="\t")

ylim <- c(min(as.numeric(Data$log2_ratio)), max(as.numeric(Data$log2_ratio)))
ylim[2] <- ylim[2]+0.5
pdf(paste(opt$prefix,".seg_plot.pdf", sep=""), height=5, width=18)
plot(as.numeric(Data$log2_ratio), pch=20, xlab='Position', ylab="Copy number", xaxt='n', ylim=ylim)
abline(v=cumsum(rle(Data$Chr)$lengths), col="red", lty=3)

if (!is.null(opt$centromereFile)) {
    cen <- read.table(opt$centromereFile, sep = '\t')
    for (j in unique(cen[,1])) {
        pos <- cen[which(cen[,1]==j)[1],3]
        index <- which(Data$Chromosome==j & Data$Start > pos)[1]
        if (!is.na(index)) {
            abline(v=index, col="darkgrey", lty=3)
        }
        text(cumsum(rle(Data$Chromosome)$lengths)-((rle(Data$Chromosome)$lengths)/2), ylim[2]-0.25)
    }
}
dev.off()



