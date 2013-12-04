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

if (!is.null(opt$centromereFile)) {
    cen <- read.table(opt$centromereFile, sep = '\t')
}

cn <- read.table(cnFile, header=T, as.is=T)
rm <- which(cn[,1] %in% c("MT", "Y"))
if (length(rm) > 0) { cn <- cn[-rm,]}
cn[which(cn[,1]=="X"),1] <- 23
cn[,1] <- as.numeric(cn[,1])
cn <- cn[order(cn[,1], cn[,2]),]
cn <- cbind(paste(cn[,1], cn[,2], sep="_"), cn[,c(1:3,7)])

cgh <- make_cghRaw(cn)
normalized <- normalize(cgh)
segmented <- segmentData(normalized, relSDlong=2, undo.splits="sdundo", undo.SD=1.5)
calls <- CGHcall(segmented, nclass=3)
excalls <- ExpandCGHcall(calls, segmented)

towrite <- cbind(fData(excalls), copynumber(excalls), segmented(excalls), calls(excalls))
colnames(towrite)[5:6] <- c("Segmented", "Calls")

write.table(towrite, file = paste(opt$prefix, ".varscan2copynumber.txt", sep=""), col.names=NA, quote=F, sep="\t")


colours <- towrite[,6]
colours[which(colours==0)] <- "black"
colours[which(colours== -1)] <- "darkred"
colours[which(colours==1)] <- "darkgreen"
colours[which(colours==2)] <- "green"
colours[which(colours==-2)] <- "red"

ylim <- c(min(as.numeric(towrite[,4])), max(as.numeric(towrite[,4])))
ylim[2] <- ylim[2]+0.5
pdf(paste(opt$prefix,".CBS.pdf", sep=""), height=5, width=18)
plot(as.numeric(towrite[,4]), pch=20, xlab='n', ylab="Copy number", xaxt='n',col=colours, ylim=ylim)
abline(v=cumsum(rle(towrite$Chr)$lengths), col="red", lty=3)

if (opt$centromereFile) {
    for (j in unique(cen[,1])) {
        pos <- cen[which(cen[,1]==j)[1],3]
        index <- towrite[which(towrite$Chromosome==j & towrite$Start > pos)[1],1]
        if (!is.na(index)) {
            abline(v=index, col="grey", lty=3)
        }
        text(cumsum(rle(towrite$Chromosome)$lengths)-((rle(towrite$Chromosome)$lengths)/2), ylim[2]-0.25)
    }
}
dev.off()

