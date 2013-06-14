#!/usr/bin/env Rscript
# cbinds together tab-delimited tables with rownames and outputs to STDOUT

arguments <- commandArgs(T);

tables <- list();
for (arg in arguments) {
    tables[[arg]] <- read.table(file = arg, sep = '\t', header = T, as.is = T, row.names = 1);
}

table.merged <- do.call(cbind, tables);

write.table(table.merged, sep = '\t', row.names = T, quote = F)
