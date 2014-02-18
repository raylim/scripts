#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));


options(warn = -1, error = quote({ traceback(); q('no', status = 1) }));

optList <- list(
                make_option("--geneCol1", default = NULL, help = "gene column 1"),
                make_option("--geneCol2", default = NULL, help = "gene column 2"),
                make_option("--sampleCol", default = NULL, help = "sample column"),
                make_option("--outDir", default = NULL, help = "Output dir"));

parser <- OptionParser(usage = "%prog [options] fusions_table", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;


if (length(arguments$args) != 1) {
    cat("Need input table\n");
    print_help(parser);
    stop();
} else if (is.null(opt$geneCol1) || is.null(opt$geneCol2) || is.null(opt$sampleCol)) {
    cat("Need column name\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outDir)) {
    cat("Need output dir\n");
    print_help(parser);
    stop();
}

fn <- arguments$args[1];
Data <- read.table(fn, sep = '\t', header = T, quote = '', comment.char = '', stringsAsFactors = F)

sData <- Data[,c(opt$sampleCol, opt$geneCol1, opt$geneCol2)]
sData <- sData[!duplicated(sData), ]

genes <- list()
for (i in 1:nrow(sData)) {
    genes[[sData[i,1]]] <- append(genes[[sData[i, 1]]], c(sData[i, 2], sData[i,3]))
}
allGenes <- unique(unlist(genes))

geneGeneM <- matrix(0, nrow = length(allGenes), ncol = length(allGenes), dimnames = list(allGenes, allGenes))
for (i in 1:nrow(sData)) {
    gene1 <- sData[i,2]
    gene2 <- sData[i,3]
    geneGeneM[gene1, gene2] <- geneGeneM[gene1, gene2] + 1
    geneGeneM[gene2, gene1] <- geneGeneM[gene2, gene1] + 1
}

x <- rowSums(geneGeneM) > 1
if (sum(x) > 1) {
    recurGeneGeneM <- geneGeneM[x, x]
    oo <- order(rowSums(geneGeneM), decreasing = T)
    recurGeneGeneM <- geneGeneM[oo, oo]
    fn <- paste(opt$outDir, "/recurGeneGene.txt", sep = "")
    write.table(recurGeneGeneM, file = fn, sep = '\t', quote = F, col.names = NA, row.names = T)
}

geneM <- sapply(genes, function (x) allGenes %in% x)
rownames(geneM) <- allGenes

if (sum(rowSums(geneM) > 1) > 1) {
    recurGeneM <- geneM[rowSums(geneM) > 1, ]
    oo <- order(rowSums(recurGeneM), decreasing = T)
    recurGeneM <- recurGeneM[oo, ]
    recurGeneM <- ifelse(recurGeneM, 1, 0)
    recurGeneM <- transform(recurGeneM, Sum = rowSums(recurGeneM))
    fn <- paste(opt$outDir, "/recurGenes.txt", sep = "")
    write.table(recurGeneM, file = fn, sep = '\t', quote = F, col.names = NA, row.names = T)
}

genePairs <- list()
for (i in 1:nrow(sData)) {
    genePairs[[sData[i,1]]] <- append(genePairs[[sData[i, 1]]], paste(sort(c(sData[i, 2], sData[i,3])), collapse = '|'))
}
allGenePairs <- unique(unlist(genePairs))
genePairM <- sapply(genePairs, function (x) allGenePairs %in% x)
rownames(genePairM) <- allGenePairs

if (sum(rowSums(genePairM) > 1) > 1) {
    recurGenePairM <- genePairM[rowSums(genePairM) > 1, ]
    oo <- order(rowSums(recurGenePairM), decreasing = T)
    recurGenePairM <- recurGenePairM[oo, ]
    recurGenePairM <- ifelse(recurGenePairM, 1, 0)
    recurGenePairM <- transform(recurGenePairM, Sum = rowSums(recurGenePairM))
    fn <- paste(opt$outDir, "/recurGenePairs.txt", sep = "")
    write.table(recurGenePairM, file = fn, sep = '\t', quote = F, col.names = NA, row.names = T)
}


