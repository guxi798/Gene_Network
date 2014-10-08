###########################################################
###	@Author		Xi Gu				###
###	@Time		July 28, 2014			###
###	@Function	Top-level Control		###
###########################################################

###########################################################################################
###  This script is a toplevel script which controls all script sections and different 	###
###  parameter versions. It requires some parameters as input, writes shell script to   ###
###  specify the jobs and submits to computer clusters. Hence, the pipeline must be run	###
###  on a computer cluster with Linux system rather than a PC.				###
###											###
###  Sample script: 1) Rscript NET.ControlFlow.Newtork.R one TRUE Pearson Signed Signed ###
###		       EM auto	soybean 0.1 1 3						###
###		    2) Rscript NET.ControlFlow.Newtork.R default soybean 0.1 1 3	###
###		    3) Rscript NET.ControlFlow.Newtork.R all soybean 0.1 1 3		###
###											###
###  Input: Some parameters, can be one of the three situations below:			###
###	    1) Multiple parameters: one. This means only one version (specified by user)###
###	       will be performed. Details of further parameters are listed below:	###
###		a) Whether log-transform data: TRUE or FALSE (default)			###
###		b) Correlation method: Pearson or Spearman (default)			###
###		c) Adjacency matrix: Signed, Unsigned (default) or Hybridy		###
###		d) TOM matrix: Signed, Unsigned (default) or FALSE			###
###		e) Cluster method: Hierarchical, Kmeans (default) or EM			###
###		f) Cluster number detection mode: auto or integer for cluster no.	###
###		g) Species name: should be consistent with the directory name		###
###		h) Coefficient Variance cutoff: value between 0 and 1			###
###		i) Expression cutoff: Value between min and max expression		###
###		m) Sample number cutoff: Value between 0 and number of samples		###
###	    2) One parameter: default. This means only one version (default) will be	###
###	       performed. Further default parameters are indicated as above.		###
###	    3) One parameter: all. This means all combinations of possible values of	###
###	       the further parameters will be tried. This will not be time-consuming,	###
###	       because all jobs will be run parallelly, unless there are not enough	###
###	       computer nodes. Hence, this mode will require a lot of computer sources.	###
###											###
###  Output: 1) Output will be generated in child scripts.				###
###########################################################################################

###########################################################################################
###  Clean and release the memory before start R					###
###########################################################################################

rm(list=ls())
gc()

###########################################################################################
###  Read in arguments, 1: runMode, 2: log-transformation, 3: correlation method,	###
###  4: adjacency matix, 5: topological overlap matrix, 6: cluster method.		###
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
runMode = tolower(args[1])

##### runMode set as one #####
if(runMode == "one"){
	log.method = tolower(args[2])
	cor.method = tolower(args[3])
	adj.method = tolower(args[4])
	TOM.method = tolower(args[5])
	cls.method = tolower(args[6])
	cls.mode = tolower(args[7])
	species = args[8]
	cut.cv = as.numeric(args[9])
	cut.exp = as.numeric(args[10])
	cut.sam = as.integer(args[11])

	cat("..........  The runMode has been set to One Version ..........\n", sep="")
	cat("..........  The Method will be ..........\n", sep="")
	cat("Log-transformation:\t\t", log.method, "\n", sep="")
	cat("Correlation Method:\t\t", cor.method, "\n", sep="")
	cat("Adjacency Matrix:\t\t", adj.method, "\n", sep="")
	cat("Topological Overlap Matrix:\t\t", TOM.method, "\n", sep="")
	cat("Clustering Method:\t\t", cls.method, "\n\n", sep="")

##### runMode set as all #####
}else if(runMode == "all"){
	log.method = c("true", "false")
	cor.method = c("pearson", "spearman")
	adj.method = c("signed", "unsigned", "hybrid")
	TOM.method = c("signed", "unsigned", "false")
	cls.method = c("hierarchical", "kmeans", "expectation-maximization")
	cls.mode = "auto"
	species = args[2]
	cut.cv = as.numeric(args[3])
	cut.exp = as.numeric(args[4])
	cut.sam = as.integer(args[5])

	cat("..........  The runMode has been set to All Versions ..........\n", sep="")
	cat("..........  All Combinations of Methods will be ..........\n", sep="")
	cat("Log-transformation:\t\t", paste(log.method, sep=", "), "\n", sep="")
	cat("Correlation Method:\t\t", paste(cor.method, sep=", "), "\n", sep="")
	cat("Adjacency Matrix:\t\t", paste(adj.method, sep=", "), "\n", sep="")
	cat("Topological Overlap Matrix:\t\t", paste(TOM.method, sep=", "), "\n", sep="")
	cat("Clustering Method:\t\t", paste(cls.method, sep=", "), "\n\n", sep="")

##### runMode set as default #####
}else if(runMode == "default"){
	log.method = "false"
	cor.method = "spearman"
	adj.method = "unsigned"
	TOM.method = "unsigned"
	cls.method = "kmeans"
	cls.mode = "auto"
	species = args[2]
	cut.cv = as.numeric(args[3])
	cut.exp = as.numeric(args[4])
	cut.sam = as.integer(args[5])

	cat("..........  The runMode has been set to Default Version ..........\n", sep="")
	cat("..........  The Method will be ..........\n", sep="")
	cat("Log-transformation:\t\t", log.method, "\n", sep="")
	cat("Correlation Method:\t\t", cor.method, "\n", sep="")
	cat("Adjacency Matrix:\t\t", adj.method, "\n", sep="")
	cat("Topological Overlap Matrix:\t\t", TOM.method, "\n", sep="")
	cat("Clustering Method:\t\t", cls.method, "\n\n", sep="")

##### runMode parameter is not acceptable #####
}else{
	stop("ERROR: the runMode has to be set one of these: ONE, ALL, DEFAULT\n\n")
}

###########################################################################################
###  Set directories for data, RData objects, result files and figures. Sub-directories	###
###  are created for RData object and Figures.						###
###########################################################################################

DataDir = "00.data"
RDataDir = "02.RData"
ResultDir = "03.result"
FigureDir = "05.graph"

##### Subdirectories under RDataDir #####
rdataNew = c("01.datExpr", "02.similarityMatrix", "03.clusterDetection", "04.moduleConversion", "05.moduleMatch", "06.networkProperties", "07.GOenrichment", "08.GSEA", "09.annotAssessment")
##### Subdirectories under RDataDir #####
resultNew = c("01.datExpr", "02.similarityMatrix", "03.clusterDetection", "04.moduleConversion", "05.moduleMatch", "06.networkProperties", "07.GOenrichment", "08.GSEA", "09.annotAssessment")
##### Subdirectories under FigureDir ####
figureNew = c("01.datExpr", "02.similarityMatrix", "03.clusterDetection", "04.moduleConversion", "05.moduleMatch", "06.networkProperties", "07.GOenrichment", "08.GSEA", "09.annotAssessment")

##### Create subdirectories ####
cat("..........  Created these directories: ..........\n", sep="")
for(i in 1:length(rdataNew)){				## for RData object
	dir.create(paste(".", species, RDataDir, rdataNew[i], sep="/"), showWarnings=FALSE, recursive=TRUE)
	cat(paste(".", species, RDataDir, rdataNew[i], sep="/"), sep="\n")
}
for(i in 1:length(resultNew)){				## for RData object
	dir.create(paste(".", species, ResultDir, resultNew[i], sep="/"), showWarnings=FALSE, recursive=TRUE)
	cat(paste(".", species, ResultDir, resultNew[i], sep="/"), sep="\n")
}
for(i in 1:length(figureNew)){				## for figures
	dir.create(paste(".", species, FigureDir, figureNew[i], sep="/"), showWarnings=FALSE, recursive=TRUE)
	cat(paste(".", species, FigureDir, figureNew[i], sep="/"), sep="\n")
}

###########################################################################################
###  Contruct functions necessary for future call.					###
###########################################################################################

buildScript <- function(log.method, 
			cor.method, 
			adj.method, 
			TOM.method, 
			cls.method,
			cls.mode, 
			species,
			cut.cv,
			cut.exp,
			cut.sam,
			queue="rcc-30d")
{
	method.name = paste(species, log.method, cor.method, adj.method, TOM.method, cls.method, cls.mode, sep=".")
	scriptFile = paste("./01.script/networkScript", method.name, "sh", sep=".")
	write("#!/bin/bash", file=scriptFile, append=FALSE, sep="\n\n")
	write(paste("time /usr/local/R/3.0.3/bin/Rscript ./01.script/NET.ConstructNetwork.R", 
		    log.method, 
		    cor.method, 
		    adj.method, 
		    TOM.method, 
		    species,
		    cut.cv,
		    cut.exp,
		    cut.sam,
		    sep=" "),
	      file = scriptFile, append = TRUE, sep = "\n")
	write(paste("time /usr/local/R/3.0.3/bin/Rscript ./01.script/NET.DetectCluster.R", 
	            log.method, 
		    cor.method, 
		    adj.method, 
		    TOM.method, 
		    cls.method, 
		    cls.mode, 
		    species,
		    sep=" "),
	      file = scriptFile, append = TRUE, sep = "\n")
	cat("\n..........  Submitted Jobs: ..........\n", sep="")
	cat(paste("qsub -q ", queue, " -e ./01.script/e.networkScript", method.name, " -o ./01.script/o.networkScript", method.name, " ", scriptFile, sep=""), sep="\n\n")
	system(paste("chmod u+x", scriptFile, sep=" "))
	system(paste("qsub -q ", queue, " -e ./01.script/e.networkScript", method.name, " -o ./01.script/o.networkScript", method.name, " ", scriptFile, sep=""))
}

cartProduct <- function(currentMatrix, newElement){
	if(length(dim(NewElement)) != 0 ){
		warning("New vector has more than one dimension.")
		return (NULL)
	}
 
	if (length(dim(currentMatrix)) == 0){
		currentRows = length(currentMatrix)
		currentMatrix = as.matrix(currentMatrix, nrow = currentRows, ncol = 1)
	} else {
		currentRows = nrow(currentMatrix)
	}
 
	var1 = replicate(length(newElement), currentMatrix, simplify=F)
	var1 = do.call("rbind", var1)
 
	var2 = rep(newElement, currentRows)
	var2 = matrix(var2[order(var2)], nrow = length(var2), ncol = 1)
 
	cartProduct = cbind(var1, var2)
	return (cartProduct)
}

###########################################################################################
###  Build script(s) and submit job(s) to the computer queue(s).			###
###########################################################################################

if(runMode == "one" || runMode == "default"){
	buildScript(log.method, 
		    cor.method, 
		    adj.method, 
		    TOM.method, 
		    cls.method, 
		    cls.mode, 
		    species,
		    cut.cv,
		    cut.exp,
		    cut.sam,
		    queue="rcc-m128-30d")
}else if(runMode == "all"){
	combo = cartProduct(log.method, cor.method)
	combo = cartProduct(combo, adj.method)
	combo = cartProduct(combo, TOM.method)
	combo = cartProduct(combo, cls.method)
	combo = cartProduct(combo, cls.mode)
	combo = cartProduct(combo, species)
	combo = cartProduct(combo, cut.cv)
	combo = cartProduct(combo, cut.exp)
	combo = cartProduct(combo, cut.sam)
	mapply(buildScript, 
		combo[,1], 
		combo[,2], 
		combo[,3], 
		combo[,4], 
		combo[,5], 
		combo[,6], 
		combo[,7],
		combo[,8],
		combo[,9],
		combo[,10],
		SIMPLIFY = FALSE)
}

###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list = ls())
gc()
