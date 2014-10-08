###########################################################
###	@Author		Xi Gu				###
###	@Time		July 24, 2014			###
###	@Function	Construct Coexpression Net	###
###########################################################

###########################################################################################
###  This script constructs gene coexpression network based on correlation coefficient,	###
###  adjacency matrix and topological overlap matrix. It takes RNAseq gene expression	###
###  file generated from cufflinks as input, filters non-expressed genes (namely zero	###
###  values in all tissues/samples) and/or low-quality genes (low values or constant	###
###  values in most tissues/samples). The criteria are taken as input parameters.	###
###  Along with the output, it also generates some statistics figures.			###
###											###
###  Input: 1) gene expression file from cuffdiff.					###
###	    2) tissue info file.							###
###	    3) detail parameters from the control script.				###
###											###
###  Output: 1) data expression in RData object, with non-expressed genes filterd.	###
###	     2) data expression in RData object, with low-quality genes filterd.	###
###	     3) topological overlap matrix in RData object if TOM is not set to FALSE.	###
###########################################################################################

###########################################################################################
###  Clean and release the memory before start R					###
###########################################################################################

rm(list=ls())
gc()

###########################################################################################
###  Load necessary R libraries and external functions for later usage.			###
###########################################################################################

cat("\n\n########## Loading Libraries ##########\n", sep="")
if(!require("gplots")){
	install.packages("gplots", dependencies = TRUE)
	library(gplots)
}
if(!require("RColorBrewer")){
	install.packages("RColorBrewer", dependencies = TRUE)
	library(RColorBrewer)
}
if(!require("devEMF")){
	install.packages("devEMF", dependencies = TRUE)
	library(devEMF)
}
if(!require("WGCNA")){
	install.packages("WGCNA", dependencies = TRUE)
	library(WGCNA)
}

source("./01.script/myfunctions.R")

###########################################################################################
###  Read in arguments.									###
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
fpkmFile = "genes.fpkm_tracking"			## gene expression file from cuffdiff
tissueFile = "record.txt"				## tissue info file

log.method = tolower(args[1])				## log-tranform: true or false
cor.method = tolower(args[2])				## correlation method: pearson or spearman
adj.method = tolower(args[3])				## adjacency method: signed, unsigned, hybrid
TOM.method = tolower(args[4])				## tom method: signed, unsigned, false
species = tolower(args[5])

cut.cv = as.numeric(args[6])				## coefficient of variance cutoff
cut.exp = as.numeric(args[7])				## expression cutoff
cut.sam = as.integer(args[8])				## at least this no of tissues have higher than cut.exp

DataDir = "00.data"
RDataDir = "02.RData"
ResultDir = "03.result"
FigureDir = "05.graph"

###########################################################################################
###  Read in fpkm using functions built in "myfunctions.R" file. cv.graph sets the 	###
###  format to generate coefficient variance (cv) distribution graph. It can take eps,	###
###  pdf, jpeg, bmp, png, tiff and emf. If this parameter is set to FALSE, no graph	###
###  will be produced where the program will be faster.					###
###########################################################################################

fpkmFile = paste(".", species, DataDir, fpkmFile, sep="/")
tissueFile = paste(".", species, DataDir, tissueFile, sep="/")
fpkm = read.fpkm(fpkmFile, tissueFile, cv.graph=FALSE)
fpkm = filter.zero(fpkm)

if(log.method == "true"){
	fpkm[fpkm == 0] = min(fpkm)
	fpkm = log(fpkm)
}

datExpr0 = filter.low(fpkm, 
		      cut.cv = 0, 
		      cut.exp = 0, 
		      cut.sam = 0)
file = paste(".", species, RDataDir, "01.datExpr", paste("01.datExpr", log.method, "QC0.RData", sep="."), sep="/")
save(datExpr0, file = file)

datExpr = filter.low(fpkm, 
		     cut.cv = cut.cv, 
		     cut.exp = cut.exp, 
		     cut.sam = cut.sam)
file = paste(".", species, RDataDir, "01.datExpr", paste("01.datExpr", log.method, "QC.RData", sep="."), sep="/")
save(datExpr, file = file)

###########################################################################################
###  Plot heatmap for all tissues/samples. Each cell value represents the correlation	###
###  between a pairwise sample pair in terms all retained genes past the filtering.	###
###########################################################################################

my.palette = colorRampPalette(c("blue", "ivory", "red"))(n=299)
figureFile = paste(".", species, FigureDir, "01.datExpr", paste("01.tissueHeatmap", log.method, "emf", sep="."), sep="/")
emf(file = figureFile, 
	width = 6, 
	height = 6, 
	family="Helvetica"
)

heatmap.2(cor(t(datExpr),method="spearman"),
	  scale = "none",
	  main = "",
	  density.info = "none", 
	  trace = "none",
	  keysize = 1,
	  margins = c(11,11),
	  col = my.palette,
	  cexRow = 1.5,
	  cexCol = 1.5
	)
dev.off()

###########################################################################################
###  Construct network using adjacency matrix and TOM.					###
###########################################################################################

###### set parameters ######
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

##### plot soft-thresholding and pick an optimal thresholding value #####
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
file = paste(".", species, RDataDir, "02.similarityMatrix", paste("02.softThresholding", log.method, "RData", sep="."), sep="/")
save(sft, file = file)
cat(paste("**Estimated soft-thresholding power: ", sft$powerEstimate), sep="\n\n")

figureFile = paste(".", species, FigureDir, "02.similarityMatrix", paste("02.softThresholding", log.method, "emf", sep="."), sep="/")
plot.sft(sft = sft, 
	 file = figureFile, 
	 format = "emf", 
	 width = 864, 
	 height = 480
	 )

##### construct network #####
if(TOM.method != "false"){
	softPower = sft$powerEstimate
	adjacency = adjacency(datExpr, 
				type = adj.method,
				power = softPower, 
				corOptions = paste("use='p', method='", cor.method, "'", sep="")
		      )

	similarityMatrix = TOMsimilarity(adjacency, TOMType = TOM.method)
	rownames(similarityMatrix) = colnames(datExpr)
	colnames(similarityMatrix) = colnames(datExpr)
	file = paste(".", species, RDataDir, "02.similarityMatrix", paste("03.TOM", log.method, cor.method, adj.method, TOM.method, "RData", sep="."), sep="/")
	save(similarityMatrix, file = file)
}else{
	softPower = 1
	similarityMatrix = adjacency(datExpr, 
				type = adj.method,
				power = softPower, 
				corOptions = paste("use='p', method='", cor.method, "'", sep="")
		      )
	rownames(similarityMatrix) = colnames(datExpr)
	colnames(similarityMatrix) = colnames(datExpr)
	file = paste(".", species, RDataDir, "02.similarityMatrix", paste("03.COR", log.method, cor.method, adj.method, TOM.method, "RData", sep="."), sep="/")
	save(similarityMatrix, file = file)
}

##### generate figures for evaluation #####
MEDissThres = 0.2
clust = TOM.cluster(similarityMatrix, 
		    clust.method = "average", 
		    cor.method = cor.method, 
		    MEDissThres = MEDissThres, 
		    minModuleSize = 100)

geneTree        = clust$geneTree
dynamicColors   = clust$dynamicColors
mergedColors    = clust$mergedColors
MEs             = clust$MEs
moduleLabels    = clust$moduleLabels
moduleColors    = clust$moduleColors
plotTOM         = clust$plotTOM

cat("..........  Created plots for eigengene clusters ..........\n", sep="")

figureFile = paste(".", species, FigureDir, "02.similarityMatrix", paste("03.EigenClust", log.method, cor.method, adj.method, TOM.method, sep="."), sep="/")
plot.cluster.eigengene(clust$METree, 
		       file = figureFile, 
		       format = "emf", 
		       width = 7, 
		       height = 6, 
		       MEDissThres = MEDissThres)

figureFile = paste(".", species, FigureDir, "02.similarityMatrix", paste("03.EigenMergedClust", log.method, cor.method, adj.method, TOM.method, sep="."), sep="/")
plot.cluster.eigengene(clust$mergedMETree, 
		       file = figureFile, 
		       format = "emf", 
		       width = 7, 
		       height = 6, 
		       MEDissThres = MEDissThres)

cat("..........  Created heatmap for similarity Matrix ..........\n", sep="")

figureFile = paste(".", species, FigureDir, "02.similarityMatrix", paste("04.similarityMatrixplot", log.method, cor.method, adj.method, TOM.method, sep="."), sep="/")
plot.TOM(plotTOM, 
	 geneTree, 
	 moduleColors, 
	 file = figureFile, 
	 format = "emf", 
	 width = 6, 
	 height = 6)

###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list = ls())
gc()


