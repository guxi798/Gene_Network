###########################################################
###	@Author		Xi Gu				###
###	@Time		July 25, 2014			###
###	@Function	Cluster Detection		###
###########################################################

###########################################################################################
###  This script detects clusters using K-means, based on the network constructed  	###
###  earlier. As the number of clusters is an important pre-set parameter for K-means,	###
###  the script will also try to optimize the cluster number automatically, unless user	###
###  specify the parameter as an input for the script. In the latter case, the script	###
###  will repress the cluster number detection module and use the parameter set by user	###
###  to detect clusters.								###
###											###
###  Input: 1) cluster Mode: integer (target cluster number) or "auto".			###
###	    2) detail parameters from the control script.				###
###											###
###  Output: 1) Clusters detected by K-means.						###
###	     2) Figures summarizing statistics about different cluster number.		###
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
if(!require("mclust")){
	install.packages("mclust", dependencies = TRUE)
	library(mclust)
}
if(!require("devEMF")){
	install.packages("devEMF", dependencies = TRUE)
	library(devEMF)
}

source("./01.script/myfunctions.R")

###########################################################################################
###  Read in arguments, the first argument indicate whether user wants to pre-specify	###
###  the number for cluster detection (any integer), or let the script determine cluster###
###  number by itself (FALSE).								###
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")

log.method = tolower(args[1])				## log-tranform: true or false
cor.method = tolower(args[2])				## correlation method: pearson or spearman
adj.method = tolower(args[3])				## adjacency method: signed, unsigned, hybrid
TOM.method = tolower(args[4])				## tom method: signed, unsigned, false
cls.method = tolower(args[5])				## cluster method
cls.mode = tolower(args[6])				## cluster mode
species = tolower(args[7])				## species folder

DataDir = "00.data"
RDataDir = "02.RData"
ResultDir = "03.result"
FigureDir = "05.graph"

###########################################################################################
###  Load RData objects, which should have been saved from upstream analyses.		###
###########################################################################################

cat("\n\n########## Loading RData Objects ##########\n", sep="")
dataFile = paste(".", species, RDataDir, "01.datExpr", paste("01.datExpr", log.method, "QC.RData", sep="."), sep="/")
	cat(".......... Loading datExpr Object ..........\n", sep="")
	cat("  ", dataFile, "\n", sep="")
	load(file = dataFile)

if(TOM.method == "false"){
	simFile = paste(".", species, RDataDir, "02.similarityMatrix", paste("03.COR", log.method, cor.method, adj.method, TOM.method, "RData", sep="."), sep="/")
	cat(".......... The Similarity Matrix is set to COR ..........\n", sep="")
	cat(".......... Loading COR Object ..........\n", sep="")
	cat("  ", simFile, "\n", sep="")
	load(file = simFile)
}else{
	simFile = paste(".", species, RDataDir, "02.similarityMatrix", paste("03.TOM", log.method, cor.method, adj.method, TOM.method, "RData", sep="."), sep="/")
	cat(".......... The Similarity Matrix is set to TOM ..........\n", sep="")
	cat(".......... Loading TOM Object ..........\n", sep="")
	cat("  ", simFile, "\n", sep="")
	load(file = simFile)
}

###########################################################################################
###  Determine the number of clusters using the Bayesian Information Criterion for	###
###  Expectation-Maximization (BIC-EM).	If clustMode is an integer, this part will be	###
###  skipped. If clustMode is set to FALSE, this part will be performed.		###
###########################################################################################

if(cls.mode == "auto"){
	cat(".......... Trying to pick up an optimal cluster number ..........\n", sep="")
	cat(".......... Bayesian Information Criterion for Expectation-maximization ..........\n", sep="")
	figureFile = paste(".", species, FigureDir, "03.clusterDetection", paste("04.ClusterNo", log.method, cor.method, adj.method, TOM.method, "emf", sep="."), sep="/")
	emf(file = figureFile, 
		width = 6, 
		height = 15, 
		family="Helvetica"
	)

	clust = Mclust(similarityMatrix, G = 1:15)
	file = paste(".", species, RDataDir, "03.clusterDetection", paste("04.ClusterNo", cls.method, log.method, cor.method, adj.method, TOM.method, "RData", sep="."), sep="/")
	save(clust, file = file)

	clustNo = dim(clust$z)[2]
	cat("model-based optimal number of clusters:", clustNo, "\n")
	plot(clust)

	dev.off()

	file = paste(".", species, RDataDir, "03.clusterDetection", paste("04.Cluster.EM", log.method, cor.method, adj.method, TOM.method, "RData", sep="."), sep="/")
	save(clust, file = file)
}else{
	cat(".......... The cluster no has been preset: ", cls.mode, " ..........\n", sep="")
	clustNo = as.integer(cls.mode)
}

###########################################################################################
###  Clustering using cluster method as specified with the no of clusters determined.	###
###########################################################################################

cat(".......... Performing Clustering ..........\n", sep="")

if(cls.method == "hierarchical"){
	if(TOM.method == "false"){
		hclust = hclust(1 - cor(similarityMatrix, method = cor.method), method = "complete")
	}else{
		hclust = hclust(1 - similarityMatrix, method = "complete")
	}
	clust = cutree(hclust, k = clustNo)
}else if (cls.method == "kmeans"){
	clust = clust.kmeans(similarityMatrix, 
			     centers = clustNo,
			     nstart = 10,
			     iter.max = 500)
}else if (cls.method == "EM"){
	clust = Mclust(similarityMatrix, G = clustNo)			## may improve to save the time
}

file = paste(".", species, RDataDir, "03.clusterDetection", paste("04.Cluster", cls.method, log.method, cor.method, adj.method, TOM.method, "RData", sep="."), sep="/")
save(clust, file = file)

###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list = ls())
gc()



