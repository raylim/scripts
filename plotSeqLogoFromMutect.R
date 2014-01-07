#!/usr/bin/env Rscript
# plot berry logos using mutect input

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("berrylogo"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outDir", default = NULL, type = "character", action = "store", help = "Output dir (required)"))

parser <- OptionParser(usage = "%prog [options]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need mutect results table\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outDir)) {
    cat("Need output directory\n\n");
    print_help(parser);
    stop();
} else {
    mutectFiles <- arguments$args;
}

results <- list()
for (mutectFile in mutectFiles) {
    sn <- sub('\\..*', '', mutectFile)
    sn <- sub('.*/', '', sn)

    d <- read.table(mutectFile, header=T, as.is = T, sep = '\t');
    d <- subset(d, judgement == "KEEP" & dbsnp_site == "NOVEL")

    d <- transform(d, mut = paste(ref_allele, alt_allele, sep = ">"))
    splitd <- split(d, d$mut)
    for (mut in names(splitd)) {
        context <- splitd[[mut]]$context
        #alt <- splitd[[mut]]$alt_allele[1]
        context <- sub("x", '', context)
        contextSplit <- t(sapply(context, function(x) unlist(strsplit(x, ''))))
        contextTable <- apply(contextSplit, 2, function(x) table(factor(x, levels = c("A", "T", "C", "G"))))
        pwm <- apply(contextTable, 2, function(x) x / sum(x))

        m <- sub('>', "to", mut)
        fn <- paste(opt$outDir, "/", sn, ".", m,  ".berry.png", sep = "")
        png(fn, height = 800, width = 800)
        berrylogo(pwm)
        dev.off()
        results[[mut]][[sn]] <- list(plotFn = fn, n = length(context))
    }
}

for (mut in names(results)) {
    m <- sub('>', "to", mut)
    pg <- openPage(paste(m, '.html', sep = ''), dirname = opt$outDir)
    for (sn in names(results[[mut]])) {
        hwrite(sn, pg, br = T)
        hwrite(paste("n =", results[[mut]][[sn]]$n), pg, br = T)
        hwriteImage(basename(results[[mut]][[sn]]$plotFn), pg)
    }
    closePage(pg)
}
