#!/usr/bin/env Rscript
# plot nucleotide frequency using mutect input

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("hwriter"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outDir", default = NULL, type = "character", action = "store", help = "Output directory (required)"))

parser <- OptionParser(usage = "%prog [options]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need mutect results table(s)\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outDir)) {
    cat("Need output dir\n\n");
    print_help(parser);
    stop();
} else {
    mutectFiles <- arguments$args;
}


d <- list()
for (mutectFile in mutectFiles) {
    cat("Reading", mutectFile, "\n")
    sp <- sub('.*/', '', mutectFile)
    sp <- sub('\\..*', '', sp)
    dat <- read.table(mutectFile, header=T, as.is = T, sep = '\t');
    dat <- subset(dat, judgement == "KEEP")
    if (!is.null(dat) &&  nrow(dat) > 0) {
        d[[sp]] <- dat
    }
}
dd <- do.call('rbind', d)

pg <- openPage('index.html', dirname = opt$outDir)

# mosaic plots for each conteingcy table
cols <- c("blue", "red", "green", "yellow")

conTable <- table(ref = dd$ref_allele, alt = dd$alt_allele)
fn <- paste(opt$outDir, '/all.mosaicplot.png', sep = '')
png(fn, height = 500, width = 500, type = 'cairo-png')
null <- mosaicplot(conTable, main = "All", col = cols)
dev.off()
hwriteImage(basename(fn), pg)
x <- as.matrix(conTable)
class(x) <- 'matrix'
hwrite(x, pg, br = T)

conTables <- lapply(d, function(x) table(ref = x$ref_allele, alt = x$alt_allele))
dfs <- lapply(conTables, as.data.frame)
df <- dfs[[1]][, c(1, 2)]
for (sp in names(dfs)) {
    x <- dfs[[sp]]
    colnames(x)[3] <- sp
    df <- merge(df, x, by = c("ref", "alt"))
}
hwrite(df, pg, br = T)
fn <- paste(opt$outDir, '/refvarCount.txt', sep = '')
write.table(df, file = fn, sep = '\t', quote = F)


for (sp in names(conTables)) {
    tab <- conTables[[sp]]
    fn <- paste(opt$outDir, '/', sp, '.mosaicplot.png', sep = '')
    png(fn, height = 400, width = 400, type = 'cairo-png')
    null <- mosaicplot(tab, col = cols, main = sp)
    dev.off()
    hwriteImage(basename(fn), pg, br = F)
    x <- as.matrix(tab)
    class(x) <- 'matrix'
    hwrite(x, pg, br = T)
}

closePage(pg)
