###########################################################
###	@Author		Xi Gu				###
###	@Time		June 12, 2014			###
###	@Function	Convert Clusters to Modules	###
###########################################################

###########################################################################################
###  This script converts clusters defined in upstream analysis into modules. Namely, 	###
###  it gives labels and colors to each cluster and assign genes to clusters. Labels	###
###  are used for calling modules and colors are for visualization.			###
###											###
###  Input: 1) expression data stored in RData object.					###
###	    2) gene clusters in RData object.						###
###	    3) whether plot or not.							###
###											###
###  Output: 1) gene module stored in RData object.					###
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

DataDir = "00.data"
RDataDir = "02.RData"
ResultDir = "03.result"
FigureDir = "05.graph"

###########################################################################################
###  Load RData objects, which should have been saved from upstream analyses.		###
###########################################################################################

cat("\n\n########## Loading RData Objects ##########\n", sep="")

exprFile = paste(".", species, RDataDir, "01.datExpr", "01.datExpr.QC.RData", sep="/")
cat("  ", exprFile, "\n", sep="")
load(file = exprFile)

clusterFile = paste(".", species, RDataDir, "03.clusterDetection", "04.Cluster.RData", sep="/")
cat("  ", clusterFile, "\n", sep="")
load(file = clusterFile)

###########################################################################################
###  Convert clusters into modules.							###
###########################################################################################

cat("\n\n########## Converting Clusters to Modules ##########\n", sep="")
cluster = clust$cluster

myColors = brewer.pal(max(cluster), "Set3")				## generate color sets
modColors = cluster
for(i in 1:max(cluster)){
	modColors[modColors == i] = myColors[i]				## assign colors to genes
}
modColors = labels2colors(modColors)					## use word to represent colors

module = list()
for(i in 1:max(cluster)){
	module = c(module, list(colnames(datExpr)[cluster == i]))		## assign genes to modules
}
names(module) = labels2colors(myColors[1:max(cluster)])			## color modules

modFile = paste(".", species, RDataDir, "04.moduleConversion", "04.networkModule.RData", sep="/")
save(module, modColors, file=modFile)					## save into RData

###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list=ls())
gc()

