#!/usr/bin/env Rscript
# rbinds together tab-delimited tables and outputs to STDOUT

options(error=traceback)

arguments <- commandArgs(T);

tables <- list();
for (arg in arguments) {
    tables[[arg]] <- read.table(file = arg, sep = '\t', header = T, as.is = T);
}

fields <- unique(unlist(c(sapply(tables, names))));

for (arg in arguments) {
    miss <- setdiff(fields, colnames(tables[[arg]]));
    tables[[arg]][,miss] <- NA;
    tables[[arg]] <- tables[[arg]][,fields];
}
table.merged <- do.call(rbind, tables);

write.table(table.merged, sep = '\t', row.names = F, quote = F)
