## default values
plate.name <- "Plate"
snp.fn <- "apt/quant-norm.pm-only.med-polish.expr.summary.txt"
genome.build <-  "hg19"
phased.bgl.dir <- "/home/ngk1/reference/phasedBGL"
calls.fn <- "apt/birdseed-v1.calls.txt"
clusters.fn <- "apt/birdseed-v1.snp-models.txt"
genome.build <- "hg19"
disease <- "breastcancer"
platform <- "SNP_6.0"
tumour <- ""
normal <- ""
results.dir <- ""

args <- (commandArgs(TRUE))

if (length(args) == 0) {
	print("No arguments supplied.\n")
    exit(1)
} else {
	for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
    }
}

if (tumour == "") {
	print("Need tumour sample name\n")
    exit(1)
}
if (results.dir == "") {
	print("Need results directory\n")
    exit(1)
}
if (out.file == "") {
	print("Need output file\n")
    exit(1)
}
if (normal != "") {
    useNormal = T
} else {
    useNormal = F
}

library(HAPSEG)
RunHapSeg(out.file = out.file, plate.name = plate.name, array.name = tumour, seg.fn = NULL, 
	snp.fn = snp.fn, 
	genome.build = genome.build, results.dir = results.dir,
	platform = platform, use.pop = "CEPH", impute.gt = TRUE, plot.segfit = TRUE, 
	merge.small = TRUE, merge.close = TRUE, min.seg.size = 3, normal = FALSE, out.p = 0.05, 
	seg.merge.thresh = 0.0001, phased.bgl.dir = phased.bgl, 
	drop.x = FALSE, drop.y = TRUE, calls.fn = calls.fn , 
	calibrate.data = TRUE, use.normal = useNormal,
    mn.sample = normal, 
	clusters.fn = clusters.fn, 
	snp.file.parser = AptSnpFileParser, clusters.file.parser = BirdseedClustersFileParser, 
	verbose = TRUE, adj.atten = 0)

