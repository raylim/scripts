#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("TitanCNA"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
        make_option("--outFile", default = NULL, type = "character", action = "store", help ="output file (required)"))

parser <- OptionParser(usage = "%prog [options] [sample params files]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need sample params files\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output file\n\n")
    print_help(parser);
    stop();
}

fns <- arguments$args

Data <- list()
for (fn in fns) {
    s <- sub(".*/", "", sub('\\..*', '', fn))
    i <- as.integer(sub(".*params_", "", sub("\\.txt", "", fn)))
    inlist <- strsplit(readLines(fn), ":\t")
    params <- lapply(inlist, tail, n = -1)
    names(params) <- lapply(inlist, head, n = 1)
    params <- lapply(params, function(x) as.numeric(unlist(strsplit(x, "[[:space:]]+"))))
    paramList <- list()
    paramList[["DbwIndex"]] <- params[["S_Dbw validity index (Both)"]]
    paramList[["normalContaminationEstimate"]] <- params[["Normal contamination estimate"]]
    paramList[["avgTumorPloidyEstimate"]] <- params[["Average tumour ploidy estimate"]]

    Data[[s]][[i]] <- paramList
}

dbwM <- matrix(NA, nrow = length(Data), ncol = length(Data[[1]]), dimnames = list(names(Data)))
avgTumorPloidyEstM <- matrix(NA, nrow = length(Data), ncol = length(Data[[1]]), dimnames = list(names(Data)))
normalContamEstM <- matrix(NA, nrow = length(Data), ncol = length(Data[[1]]), dimnames = list(names(Data)))
for (s in names(Data)) {
    for (i in 1:length(Data[[1]])) {
        dbwM[s, i] <- Data[[s]][[i]]$DbwIndex
        avgTumorPloidyEstM[s, i] <- Data[[s]][[i]]$avgTumorPloidyEstimate
        normalContamEstM[s, i] <- Data[[s]][[i]]$normalContaminationEstimate
    }
}
dbwMin <- apply(dbwM,1, which.min)
avgTumorPloidyEstMin <- avgTumorPloidyEstM[cbind(1:nrow(dbwM), dbwMin)]
normalContamEstMin <- normalContamEstM[cbind(1:nrow(dbwM), dbwMin)]

results <- data.frame(row.names = names(dbwMin), argminSDbw = dbwMin, avgTumorPloidyEst = avgTumorPloidyEstMin, normalContamEst = normalContamEstMin)

write.table(dbwM, file = 'titan_dbw_index.txt', sep = '\t', quote = F, col.names = F)
write.table(results, file = 'titan_param_summary.txt', sep = '\t', quote = F)
