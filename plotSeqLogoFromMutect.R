#!/usr/bin/env Rscript
# plot berry logos using mutect input

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("berrylogo"));
suppressPackageStartupMessages(library("ggplot2"));
suppressPackageStartupMessages(library("hwriter"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--lowAF", default = F, action = "store_true", help = "only output low tumor allelic frequency variants"),
                make_option("--highAF", default = F, action = "store_true", help = "only output high tumor allelic frequency variants"),
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

contexts <- list()
results <- list()
for (mutectFile in mutectFiles) {
    sn <- sub('\\..*', '', mutectFile)
    sn <- sub('.*/', '', sn)

    cat("reading", mutectFile, "\n")
    d <- read.table(mutectFile, header=T, as.is = T, sep = '\t');
    d <- subset(d, judgement == "KEEP" & dbsnp_site == "NOVEL")
    if (opt$lowAF) {
        d <- subset(d, tumor_f < 0.15)
    }
    if (opt$highAF) {
        d <- subset(d, tumor_f > 0.15)
    }


    d <- transform(d, mut = paste(ref_allele, alt_allele, sep = ">"))
    splitd <- split(d, d$mut)
    for (mut in names(splitd)) {
        if (nrow(splitd[[mut]]) < 5) {
            next();
        }
        m <- sub('>', "to", mut)

        context <- splitd[[mut]]$context
        contextSplit <- t(sapply(context, function(x) unlist(strsplit(x, ''))))

        bpfn <- paste(opt$outDir, "/", sn, "_", m,  "_barplot1.png", sep = "")
        cat("plotting", bpfn, "\n")
        png(bpfn, height = 800, width = 800, type = 'cairo-png')
        barplot(table(as.factor(apply(contextSplit[,3:5], 1, paste, collapse = ""))), horiz = T, las = 2)
        dev.off()

        bp2fn <- paste(opt$outDir, "/", sn, "_", m,  "_barplot2.png", sep = "")
        cat("plotting", bp2fn, "\n")
        png(bp2fn, height = 800, width = 800, type = 'cairo-png')
        barplot(table(as.factor(apply(contextSplit[,3:6], 1, paste, collapse = ""))), horiz = T, las = 2)
        dev.off()
        
        bp3fn <- paste(opt$outDir, "/", sn, "_", m,  "_barplot3.png", sep = "")
        cat("plotting", bp3fn, "\n")
        png(bp3fn, height = 800, width = 800, type = 'cairo-png')
        barplot(table(as.factor(apply(contextSplit[,3:7], 1, paste, collapse = ""))), horiz = T, las = 2)
        dev.off()

        bp4fn <- paste(opt$outDir, "/", sn, "_", m,  "_barplot4.png", sep = "")
        cat("plotting", bp4fn, "\n")
        png(bp4fn, height = 800, width = 800, type = 'cairo-png')
        barplot(table(as.factor(apply(contextSplit[,c(3,4,5,7)], 1, paste, collapse = ""))), horiz = T, las = 2)
        dev.off()

        context <- sub("x", '', context)
        contexts[[mut]] <- c(contexts[[mut]], context)
        contextSplit <- t(sapply(context, function(x) unlist(strsplit(x, ''))))
        contextTable <- apply(contextSplit, 2, function(x) table(factor(x, levels = c("A", "T", "C", "G"))))
        pwm <- apply(contextTable, 2, function(x) x / sum(x))

        fn <- paste(opt$outDir, "/", sn, "_", m,  "_berry.png", sep = "")
        cat("plotting", fn, "\n")
        png(fn, height = 800, width = 800, type = 'cairo-png')
        print(berrylogo(pwm))
        dev.off()
        results[[mut]][[sn]] <- list(plotFn = fn, bpFn1 = bpfn, bpFn2 = bp2fn, bpFn3 = bp3fn, bpFn4 = bp4fn, n = length(context))
    }
}

for (mut in names(results)) {
    m <- sub('>', "to", mut)
    pg <- openPage(paste(m, '.html', sep = ''), dirname = opt$outDir)

    # plot logo using all
    contextSplit <- t(sapply(contexts[[mut]], function(x) unlist(strsplit(x, ''))))
    contextTable <- apply(contextSplit, 2, function(x) table(factor(x, levels = c("A", "T", "C", "G"))))
    pwm <- apply(contextTable, 2, function(x) x / sum(x))
    m <- sub('>', "to", mut)
    fn <- paste(opt$outDir, "/", m, "_berry.png", sep = "")
    cat("plotting", fn, "\n")
    png(fn, height = 800, width = 800, type = 'cairo-png')
    print(berrylogo(pwm))
    dev.off()
    img <- hwriteImage(basename(fn))
    hwrite(c(img, caption = "All"), pg, br = T, dim = c(2, 1), row.style = list(caption='text-align:center;background-color:#fac'), row.names = F)

    for (sn in names(results[[mut]])) {
        capt <- paste(sn, " logo (n = ", results[[mut]][[sn]]$n, ")", sep = "")
        img <- hwriteImage(basename(results[[mut]][[sn]]$plotFn))
        hwrite(c(img, caption = capt), pg, br = T, dim = c(2, 1), row.style = list(caption='text-align:center;background-color:#fac'), row.names = F)

        capt <- paste(sn, " -1 to 1 barplot (n = ", results[[mut]][[sn]]$n, ")", sep = "")
        img <- hwriteImage(basename(results[[mut]][[sn]]$bpFn1))
        hwrite(c(img, caption = capt), pg, br = T, dim = c(2, 1), row.style = list(caption='text-align:center;background-color:#fac'), row.names = F)

        capt <- paste(sn, " -1 to 2 barplot (n = ", results[[mut]][[sn]]$n, ")", sep = "")
        img <- hwriteImage(basename(results[[mut]][[sn]]$bpFn2))
        hwrite(c(img, caption = capt), pg, br = T, dim = c(2, 1), row.style = list(caption='text-align:center;background-color:#fac'), row.names = F)

        capt <- paste(sn, " -1 to 3 barplot (n = ", results[[mut]][[sn]]$n, ")", sep = "")
        img <- hwriteImage(basename(results[[mut]][[sn]]$bpFn3))
        hwrite(c(img, caption = capt), pg, br = T, dim = c(2, 1), row.style = list(caption='text-align:center;background-color:#fac'), row.names = F)

        capt <- paste(sn, " -1,1,3 barplot (n = ", results[[mut]][[sn]]$n, ")", sep = "")
        img <- hwriteImage(basename(results[[mut]][[sn]]$bpFn4))
        hwrite(c(img, caption = capt), pg, br = T, dim = c(2, 1), row.style = list(caption='text-align:center;background-color:#fac'), row.names = F)
    }
    closePage(pg)
}
