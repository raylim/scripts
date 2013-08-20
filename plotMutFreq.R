#!/usr/bin/env Rscript
# plot nucleotide frequency using mutect input

suppressPackageStartupMessages(library("optparse"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--out", default = NULL, type = "character", action = "store", help = "Output prefix (required)"))

parser <- OptionParser(usage = "%prog [options]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need mutect results table\n");
    print_help(parser);
    stop();
} else if (is.null(opt$out)) {
    cat("Need output file\n\n");
    print_help(parser);
    stop();
} else {
    mutectFiles <- arguments$args;
}

d <- list()
for (mutectFile in mutectFiles) {
    cat("Reading", mutectFile, "\n")
    d[[mutectFile]] <- read.table(mutectFile, header=T, as.is = T, sep = '\t');
    d[[mutectFile]] <- subset(d[[mutectFile]], judgement == "KEEP" & dbsnp_site == "NOVEL")
}
d <- do.call('rbind', d)


bases <- c("A", "T", "C", "G")
# ref > x
for (base in bases) {
    subd <- subset(d, ref_allele == base)
    context <- t(sapply(subd$context, function(x) unlist(strsplit(x, ''))))
    context[,4] <- base
    contextTable <- apply(context, 2, function(x) table(factor(x, levels = bases)))

    fn <- paste(opt$out, ".", base, "toX.base_freq.png", sep = "")
    png(fn, height = 800, width = 800)
    cols <- c("green", "blue", "grey", "red")
    barplot(contextTable, col = cols, space = 0, ylab = "Number of Bases", names.arg = -3:3, xlab = "Position")
    legend('top', rownames(contextTable), fill = cols)
    dev.off()
}
