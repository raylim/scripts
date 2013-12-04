#!/usr/bin/env Rscript
# segment copy number data and generate plot

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("DNAcopy"));
suppressPackageStartupMessages(library("SMAP"));
suppressPackageStartupMessages(library("CGHcall"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--prefix", default = NULL, type = "character", action = "store", help ="Output prefix (required)"))

parser <- OptionParser(usage = "%prog [options] copynumFile", option_list = optList);
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

cat("Reading", cnFile, " ... ")
cn <- read.table(cnFile, header = T)
cat("done\n")
CNA.object <- CNA(genomdat = cn$adjusted_log_ratio, chrom = cn$chrom, maploc = cn$chr_start, data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
cat("Segmenting ... ")
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
# reorder chroms
xx <- as.character(segs$output$chrom)
xx <- sub('[A-Za-z]*', '', xx, perl = T)
xx <- as.integer(xx)
oo <- order(as.integer(as.character(xx)), na.last = T)
segs$output <- segs$output[oo, ]
cat("Done\n")

fn <- paste(opt$prefix, '.segment.txt', sep = '')
cat("Writing to", fn, " ... ")
write.table(segs$output[ ,-1], file = fn , row.names=F, quote=F, sep="\t")
cat("done\n")

fn <- paste(opt$prefix, '.cnv_full.png', sep = '')
cat("Plotting to", fn, " ... ")
png(filename = fn, res = 70, width = 1000,
    height = 800, pointsize = 16, type = "cairo-png")
plot(segs)
null <- dev.off()
cat("done\n")

fn <- paste(opt$prefix, '.cnv_by_chr.png', sep = '')
cat("Plotting to", fn, " ... ")
png(filename = fn, res = 70, width = 2000,
    height = 2400, pointsize = 16, type = "cairo-png")
plot(segs, plot.type = "samplebychrom")
null <- dev.off()
cat("done\n")


# reorder chroms
xx <- as.character(cn$chrom)
xx <- sub('[A-Za-z]*', '', xx, perl = T)
xx <- as.integer(xx)
oo <- order(as.integer(as.character(xx)), na.last = T)
cn <- cn[oo, ]
obs <- SMAPObservations(value = as.numeric(cn$adjusted_log_ratio), chromosome = as.character(cn$chrom), startPosition = as.numeric(cn$chr_start), endPosition = as.numeric(cn$chr_stop))
init.means <- c(0, -1, 1, -2, 2)
init.sds <- rep(0.1, 5)
phi <- cbind(init.means, init.sds)
hmm <- SMAPHMM(noStates=5, Phi=phi, initTrans=0.02)
profile <- smap(hmm, obs, verbose=2)

fn <- paste(opt$prefix, '.cnv_full_smap.png', sep = '')
cat("Plotting to", fn, " ... ")
png(filename = fn, res = 70, width = 2000,
    height = 800, pointsize = 16, type = "cairo-png")
plot(profile, ylab = "logratio")
null <- dev.off()
cat("done\n")

fn <- paste(opt$prefix, '.cnv_full_smap_byindex.png', sep = '')
cat("Plotting to", fn, " ... ")
png(filename = fn, res = 70, width = 2000,
    height = 800, pointsize = 16, type = "cairo-png")
chrInd <- which(!duplicated(chromosome(obs)))
chrMid <- chrInd + diff(c(chrInd, length(value(obs)))) / 2
chr <- chromosome(obs)[!duplicated(chromosome(obs))]

cols <- Q(profile)
cols[cols == 1] <- 'Black' # neutral
cols[cols == 3] <- 'Darkgreen' # gain
cols[cols == 2] <- 'Darkred' # loss
cols[cols == 5] <- 'green' # amp
cols[cols == 4] <- 'Red' # deletion

plot(value(obs), col = cols, ylab = "logratio", xaxt = 'n', xlab = '', pch = 20)
axis(1, at = chrMid, labels = chr, tick = F)
abline(v = chrInd, lty = 3, col = 'grey')
null <- dev.off()
cat("done\n")


cat("Finished\n")

