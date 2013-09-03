#!/usr/bin/env Rscript
# rbinds together tab-delimited tables and outputs to STDOUT

suppressPackageStartupMessages(library("optparse"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--sampleName", action = "store_true", default = F, help = "add samplename column [default %default]"),
                make_option("--tumorNormal", action = "store_true", default = F, help = "add tumor-normal samplename column [default %default]"))

parser <- OptionParser(usage = "%prog [options] vcf.file", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;
files <- arguments$args;

tables <- list();
headers <- list();
for (f in files) {
    X <- read.table(file = f, sep = '\t', as.is = T, comment.char = '', quote = '');
    if (nrow(X) <= 1) {
        next
    }
    tables[[f]] <- X
    headers[[f]] <- tables[[f]][1,];
    headers[[f]] <- sub('#', '', headers[[f]])
    colnames(tables[[f]]) <- headers[[f]];
    tables[[f]] <- tables[[f]][-1, ];

    if (opt$sampleName) {
        sname <- sub('\\..*', '', f)
        sname <- sub('.*/', '', sname)
        tables[[f]][,"SAMPLE"] <- sname
        headers[[f]] <- c(headers[[f]], "SAMPLE")
        headers[[f]] <- sub(paste(sname, '\\.', sep = ''), 'SAMPLE.', headers[[f]])
        colnames(tables[[f]]) <- headers[[f]];
    }
    if (opt$tumorNormal) {
        sname <- sub('\\..*', '', f)
        sname <- sub('.*/', '', sname)
        tumor <- sub('_.*', '', sname)
        normal <- sub('.*_', '', sname)
        tables[[f]][,"TUMOR_SAMPLE"] <- tumor
        headers[[f]] <- c(headers[[f]], "TUMOR_SAMPLE")
        tables[[f]][,"NORMAL_SAMPLE"] <- normal
        headers[[f]] <- c(headers[[f]], "NORMAL_SAMPLE")

        headers[[f]] <- sub(paste(tumor, '\\.', sep = ''), 'TUMOR.', headers[[f]])
        headers[[f]] <- sub(paste(normal, '\\.', sep = ''), 'NORMAL.', headers[[f]])
        colnames(tables[[f]]) <- headers[[f]];
    }
}
if (length(tables) == 0) {
    quit(save = 'no', status = 0)
}

fields <- unique(unlist(headers))

for (f in names(tables)) {
    miss <- setdiff(fields, colnames(tables[[f]]));
    tables[[f]][,miss] <- NA;
    tables[[f]] <- tables[[f]][, fields];
}
table.merged <- do.call(rbind, tables);
rownames(table.merged) <- NULL

if (opt$sampleName) {
    x <- which(colnames(table.merged) == "SAMPLE")
    y <- which(colnames(table.merged) != "SAMPLE")
    table.merged <- table.merged[, c(x,y)]
}
if (opt$tumorNormal) {
    xx <- colnames(table.merged) == "TUMOR_SAMPLE" | colnames(table.merged) == "NORMAL_SAMPLE"
    x <- which(xx)
    y <- which(!xx)
    table.merged <- table.merged[, c(x,y)]
}

write.table(table.merged, sep = '\t', row.names = F, quote = F)
