###########################################################
###	@Author		Xi Gu				###
###	@Time		June 13, 2014			###
###	@Function	Compute Eigengene		###
###########################################################

###########################################################################################
###  This script computes eigengene for each module defined in upstream analysis,	###
###  either based on expression data (datExpr) or Topological Overlap Matrix (TOM).	###
###  This part is necessary because computing eigengenes based on TOM is very		###
###  computational intensive. Running the analysis once and storing results can save	###
###  considerable amount of time.							###
###											###
###  Input: 1) expression data stored in RData object.					###
###	    2) gene modules stored in RData object.					###
###	    3) TOM matrix stored in RData object.					###
###											###
###  Output: 1) Module eigengene and ordered eigengenes stored in RData object.		###
###########################################################################################

###########################################################################################
###  Clean and release the memory before start R					###
###########################################################################################

rm(list=ls())
gc()

###########################################################################################
###  Load necessary R libraries for later usage.					###
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
if(!require("psych")){
        install.packages("psych", dependencies = TRUE)
        library(psych)
}

###########################################################################################
###  Read in arguments.									###
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")

species = tolower(args[1])				## species folder
similarMethod = args[2]

DataDir = "00.data"
RDataDir = "02.RData"
ResultDir = "03.result"
FigureDir = "05.graph"

###########################################################################################
###  Load RData objects, which should have been saved from upstream analyses.		###
###########################################################################################

cat("\n\n########## Loading RData Objects ##########\n", sep="")
dataFile = paste(".", species, RDataDir, "01.datExpr", "01.datExpr.QC.RData", sep="/")
cat("  ", dataFile, "\n", sep="")
load(file = dataFile)

simFile = paste(".", species, RDataDir, "02.similarityMatrix", "03.similarityMatrix.RData", sep="/")
cat("  ", simFile, "\n", sep="")
load(file = simFile)

modFile = paste(".", species, RDataDir, "04.moduleConversion", "04.networkModule.RData", sep="/")
cat("  ", modFile, "\n", sep="")
load(modFile)

###########################################################################################
###  Compute Eigengene for each module, either based on datExpr or TOM.			###
###########################################################################################

cat("\n\n########## Computing Eigengene for the similarityMatrix ##########\n", sep="")
cat(".......... Based on datExpr ..........\n", sep="")
MEs = moduleEigengenes(datExpr, modColors)$eigengenes
MET = orderMEs(MEs)
MEfile = paste(".", species, RDataDir, "04.moduleConversion", "04.moduleEigengene.datExpr.RData", sep="/")
save(MEs, MET, file=MEfile)

tissueMEs = MEs
tissueMET = MET

cat(".......... Based on TOM ..........\n", sep="")
MEs = moduleEigengenes(similarityMatrix, modColors)$eigengenes
MET = orderMEs(MEs)
MEfile = paste(".", species, RDataDir, "04.moduleConversion", "04.moduleEigengene.similarityMatrix.RData", sep="/")
save(MEs, MET, file=MEfile)


###########################################################################################
###  Generate plot based on similarity matrix pre-specified.				###
###########################################################################################

	cat(".......... Plotting Figures ..........\n", sep="")

	## plot ME tree and Module relationship
	figureFile = paste(".", species, FigureDir, "04.moduleConversion", "04.moduleEigene.Dendrogram.emf", sep="/")
	emf(file = figureFile, 
		width = 6, 
		height = 10, 
		family="Helvetica"
	)
	cat(".......... ", figureFile, "\n", sep="")
	par(cex = 1.0)
	plotEigengeneNetworks(MET, 
			      setLabels = "Eigengene dendrogram", 
			      marDendro = c(0,4,2,4.5), 
			      marHeatmap = c(3,3.5,2,2), 
			      plotHeatmaps = TRUE)
	dev.off()

	## plot modules versus tissues
	figureFile = paste(".", species, FigureDir, "04.moduleConversion", "04.moduleEigene.Tissues.emf", sep="/")
	emf(file = figureFile,
		width = 6,
		height = 8,
		family="Helvetica"
	)
	cat(".......... ", figureFile, "\n", sep="")
	
	rownames(tissueMET) = rownames(datExpr)
	my.palette = colorRampPalette(c("blue", "ivory", "red"))(n=299)

	heatmap.2(t(t(tissueMET)),
		Rowv = TRUE,
		Colv = FALSE,
		dendrogram = "none",
		keysize = 1,
		symkey = F,
		trace = "none",
		density.info = "none",
		scale = "none",
		col = my.palette,
		srtCol = 50,
		margins=c(8,11))
	dev.off()

###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list=ls())
gc()
