#!/usr/bin/env Rscript
# compute cross correlation of two expression matrices

suppressPackageStartupMessages(library(optparse));
suppressPackageStartupMessages(library(qvalue));

usage <- "%prog [options] matrix1.file matrix2.file";

optList <- list(
                make_option("--parallelBackend", default = "mpi", help = "Parallel backend: mpi or mc  [default %default]"),
                make_option("--outputDir", default = "results", help = "Output directory [default %default]"), 
                make_option("--numCores", default = 5, help = "Number of cores to use with mc backend [default %default]"));

parser <- OptionParser(usage = usage, option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 2) {
    cat(usage, "\n");
    stop();
}

m1File <- arguments$args[1];
m2File <- arguments$args[2];
#outFile <- arguments$args[3];
outDir <- opt$outputDir;

exprs1 <- as.matrix(read.table(m1File, sep = '\t', header = T, as.is = T, row.names = 1));
exprs2 <- as.matrix(read.table(m2File, sep = '\t', header = T, as.is = T, row.names = 1));

if (opt$parallelBackend == "mpi") {
    suppressPackageStartupMessages(library(doMPI));
    cl <- startMPIcluster();
    registerDoMPI(cl);
    cat("Using", mpi.universe.size(), "threads\n");
} else {
    suppressPackageStartupMessages(library(doMC));
    registerDoMC(opt$numCores);
}

for (i in 1:(ceiling(nrow(exprs1) / 1000))) {
    dn <- sprintf('%s/%03d', outDir, i);
    dir.create(dn);
}

chunkSize <- 1;
mpi.opts <- list(chunkSize = chunkSize);
mc.opts <- list(chunkSize = chunkSize);

#rhos <- matrix(NA, nrow = nrow(exprs1), ncol = nrow(exprs2));
#pvals <- matrix(NA, nrow = nrow(exprs1), ncol = nrow(exprs2));
files <- c();
X <- foreach(i = 1:nrow(exprs1), .options.mpi = mpi.opts, .options.mc = mc.opts) %dopar% {
        #foreach(j = 1:nrow(exprs2), .combine = 'c') %dopar% {
            #cor(exprs1[i, ], exprs2[j, ]);
    cortests <- apply(exprs2, 1, function(x) cor.test(exprs1[i, ], x, method = 'spearman'));
    p <- sapply(cortests, function(x) x$p.value);
    adj.p <- qvalue(p)$qvalues;
    r <- sapply(cortests, function(x) x$estimate);
    #rhos[i, ] <- r;
    #pvals[i, ] <- p;
    fn <- sprintf('%05d.txt', i);
    dn <- sprintf('%03d', ceiling(i / 1000));
    files <- c(files, paste(dn, fn, sep = '/'));
    res <- cbind(r, p, adj.p);
    write.table(res, file = paste(outDir, dn, fn, sep = '/'), quote = F, row.names = F, col.names = F);
    #cat('Writing to', fn, "\n");
    1
    #rhos
};

write.table(rownames(exprs1), file = paste(outDir, '/', 'genes1.txt', sep = ''), quote = F, row.names = F, col.names = F);
write.table(rownames(exprs2), file = paste(outDir, '/', 'genes2.txt', sep = ''), quote = F, row.names = F, col.names = F);
write.table(files, file = paste(outDir, '/', 'files.txt', sep = ''), quote = F, row.names = F, col.names = F);

# rhos <- foreach(i = 1:nrow(exprs1), .combine = "cbind") %:% 
# foreach (j = 1:nrow(exprs2), .combine = 'c') %dopar% {
# cor(exprs1[i, ], exprs2[j, ]);
# };

#write.table(rhos, file = outFile, quote = F, sep = '\t', row.names = rownames(exprs2), col.names = rownames(exprs1));
#write.table(pvals, file = paste(outdir, '/pairCorPvals.txt', sep = ''), quote = F, sep = '\t', row.names = rownames(exprs2), col.names = rownames(exprs1));

if (opt$parallelBackend == "mpi") {
    closeCluster(cl);
    mpi.quit();
}

