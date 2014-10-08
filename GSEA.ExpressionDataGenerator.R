###########################################################
###	@Author		Xi Gu				###
###	@Time		June 05, 2014			###
###	@Function	Reformat Exprssion File		###
###########################################################

###########################################################################################
###  This script generates expression files in txt format and phenotype file in cls	###
###  format, which are required for GSEA analysis. One file per module which should 	###
###  have been obtained from upstream analysis. Each file contains expression values 	###
###  for genes assigned to the corresponding module. One row per gene. The script also	###
###  generates phenotype file which contains values for	eigengene of the corresponding	###
###  module. The phenotype files will be used in GSEA.					###
###											###
###  Input: 1) expression data stored in RData object.					###
###	    2) gene to module assignment stored in RData object.			###
###											###
###  Output: 1) expression files for each module in txt format.				###
###	     2) phenotype files for each module in cls format.				###
###  Note, both txt and cls formats are defined by GSEA and required in GSEA analysis.	###
###  Here, cls is a continuous phenotype file.						###
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
if(!require("WGCNA")){
	install.packages("WGCNA", dependencies = TRUE)
	library(WGCNA)
}

###########################################################################################
###  Read in arguments. The only argument is to choose the similarity matrix which can	###
###  either be the direct expression data or topological overlap index (TOM).		###
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
similarMethod = args[1]
species = args[2]

RDataDir = "02.RData"
ResultDir = "03.result"
FigureDir = "05.graph"

###########################################################################################
###  Load RData objects, which should have been saved from upstream analyses.		###
###########################################################################################

cat("\n\n########## Loading RData Objects ##########\n", sep="")
modFile = paste(".", species, RDataDir, "04.moduleConversion", "04.networkModule.RData", sep="/")
cat("  ", modFile, "\n", sep="")
load(modFile)

if(similarMethod == "datExpr"){
	simFile = paste(".", species, RDataDir, "01.datExpr", "01.datExpr.QC.RData", sep="/")
	load(file = simFile)
	similarityMatrix = datExpr
}else if(similarMethod == "similarityMatrix"){
	simFile = paste(".", species, RDataDir, "02.similarityMatrix", "03.similarityMatrix.RData", sep="/")
	load(file = simFile)
}
cat("  ", simFile, "\n", sep="")

###########################################################################################
###  Output expression data into required format.					###
###########################################################################################

cat("\n\n########## Writing Expression Files ##########\n", sep="")
cat(".......... Based on Similarity Matrix: ", similarMethod, " ..........\n", sep="")

genenames = colnames(similarityMatrix)
for(i in 1:length(module)){
	loc = match(module[[i]], genenames)			## find which genes belong to the module
	subSimilarityMatrix = t(similarityMatrix[,loc])		## fish corresponding values for the gene
	data = cbind(genenames[loc], rep(NA, times=length(module[[i]])), subSimilarityMatrix)
	colnames(data)[1:2] = c("NAME", "DESCRIPTION")		## the first two columns required by GSEA txt format
	file = paste(".", species, ResultDir, "08.GSEA", paste("02.ExpressionSet" , similarMethod, names(module)[i], "txt", sep="."), sep="/")
	cat(paste(".", species, ResultDir, "08.GSEA", paste("02.ExpressionSet" , similarMethod, names(module)[i], "txt", sep="."), sep="/"), sep="\n")
	write.table(data, file, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}


###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list=ls())
gc()
