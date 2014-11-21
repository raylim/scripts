#!/usr/bin/env Rscript
# plots emu signature

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("VariantAnnotation"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("hwriter"));
suppressPackageStartupMessages(library("org.Hs.eg.db"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outDir", default = NULL, help = "output dir (required)"),
                make_option("--samples", default = NULL, help = "samples file"),
                make_option("--sampleSubset", default = NULL, help = "sample subset file"),
                make_option("--inPrefix", default = NULL, help = "input prefix (required)"))

parser <- OptionParser(usage = "%prog [options]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outDir)) {
    cat("Need html output file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$inPrefix)) {
    cat("Need input prefix\n");
    print_help(parser);
    stop();
} else {
    files <- arguments$args;
}

glob <- paste(opt$inPrefix, '_*_ml_spectra.txt', sep = '')
spectraFiles <- Sys.glob(glob)

glob <- paste(opt$inPrefix, '_*_map_activities.txt', sep = '')
activityFiles <- Sys.glob(glob)

glob <- paste(opt$inPrefix, '_*_assigned.txt', sep = '')
assignedFiles <- Sys.glob(glob)

pg <- openPage('index.html', dirname = opt$outDir, title = 'EMu results')

set.seed(001)
palette(sample(rainbow(30)))

for (fn in spectraFiles) {
    spectra <- read.table(fn, sep = ' ')
    spectra <- spectra[,-97] # remove empty col
    for (i in 1:nrow(spectra)) {
        ofn <- paste(opt$outDir, "/", basename(fn), sep = '')
        ofn <- sub('\\.txt$', paste("_", i, '.pdf', sep = ''), ofn)
        pdf(ofn, height = 5, width = 10)
        cols <- rep(c('LightBlue', 'Black', 'Red', 'Grey', 'Green', 'Magenta'), each = 16)
        barplot(t(spectra[i,]) * 100, beside = T, col = cols, border = cols, xaxt = 'n', main = paste("Signature", i), col.main = i, ylab = "% of mutations")
        labs <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
        mtext(labs, side = 1, at = 1:6 * 16 - 7.5)
        null <- dev.off()
    }

    for (i in 1:nrow(spectra)) {
        ofn <- paste(opt$outDir, "/", basename(fn), sep = '')
        ofn <- sub('\\.txt$', paste("_", i, '.png', sep = ''), ofn)
        png(ofn, height = 500, width = 800, type = 'cairo')
        cols <- rep(c('LightBlue', 'Black', 'Red', 'Grey', 'Green', 'Magenta'), each = 16)
        barplot(t(spectra[i,]) * 100, beside = T, col = cols, border = cols, xaxt = 'n', main = paste("Signature", i), ylab = "% of mutations", col.main = i)
        labs <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
        mtext(labs, side = 1, at = 1:6 * 16 - 7.5)
        null <- dev.off()
        hwriteImage(basename(ofn), pg, br = T)
    }
}


samples <- scan('mutations.txt.samples', what = 'character')
sampleSubset <- scan(opt$sampleSubset, what = 'character')
x <- which(samples %in% sampleSubset)
for (fn in assignedFiles) {
    assigned <- read.table(fn, sep = ' ')
    assigned <- as.matrix(assigned[,-ncol(assigned)])
    assigned <- assigned[x, ]
    par(mar = c(5, 10, 5, 5))
    barplot(t(assigned / rowSums(assigned)), col = 1:5, space = 0, border = F, horiz = T, names = samples[x], las = 2)
}

closePage(pg)
