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
###  Sample script: 1) Rscript NET.ControlFlow.Cluster.R one TRUE Pearson Signed Signed ###
###		       EM auto soybean							###
###		    2) Rscript NET.ControlFlow.Cluster.R default soybean		###
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
	species = unlist(strsplit(args[8], ","))

##### runMode set as default #####
}else if(runMode == "default"){
	log.method = "false"
	cor.method = "spearman"
	adj.method = "unsigned"
	TOM.method = "unsigned"
	cls.method = "kmeans"
	cls.mode = "auto"
	species = unlist(strsplit(args[2], ","))

##### runMode parameter is not acceptable #####
}else{
	stop("ERROR: the runMode has to be set one of these: ONE, DEFAULT\n\n")
}

###########################################################################################
###  Set directories for data, RData objects, result files and figures. Sub-directories	###
###  are created for RData object and Figures.						###
###########################################################################################

DataDir = "00.data"
RDataDir = "02.RData"
ResultDir = "03.result"
FigureDir = "05.graph"

###########################################################################################
###  Rename picked version to the standard version for downstream analysis.		###
###########################################################################################



###########################################################################################
###  Contruct functions necessary for future call.					###
###########################################################################################

buildScript <- function(sp, queue="rcc-30d")
{
	##### datExpr file #####
	oldFile = paste(".", sp, RDataDir, "01.datExpr", paste("01.datExpr", log.method, "QC", "RData", sep="."), sep="/")
	newFile = paste(".", sp, RDataDir, "01.datExpr", "01.datExpr.QC.RData", sep="/")
	system(paste("cp", oldFile, newFile, sep=" "))
	
	oldFile = paste(".", sp, RDataDir, "01.datExpr", paste("01.datExpr", log.method, "QC0", "RData", sep="."), sep="/")
	newFile = paste(".", sp, RDataDir, "01.datExpr", "01.datExpr.QC0.RData", sep="/")
	system(paste("cp", oldFile, newFile, sep=" "))
	
	##### similarity matrix file #####
	if(TOM.method == "false"){
		sim = "COR"
	}else{
		sim = "TOM"
	}
	
	oldFile = paste(".", sp, RDataDir, "02.similarityMatrix", 
			paste("03", sim, log.method, cor.method, adj.method, TOM.method, "RData", sep="."), 
			sep="/")
	newFile = paste(".", sp, RDataDir, "02.similarityMatrix", "03.similarityMatrix.RData", sep="/")
	system(paste("cp", oldFile, newFile, sep=" "))
	
	##### cluster file #####
	oldFile = paste(".", sp, RDataDir, "03.clusterDetection", 
			paste("04.Cluster", cls.method, log.method, cor.method, adj.method, TOM.method, "RData", sep="."), 
			sep="/")
	newFile = paste(".", sp, RDataDir, "03.clusterDetection", "04.Cluster.RData", sep="/")
	system(paste("cp", oldFile, newFile, sep=" "))
	
	##### write shell script #####
	scriptFile = paste("./01.script/clusterScript", sp, "sh", sep=".")
	write("#!/bin/bash", file=scriptFile, append=FALSE, sep="\n\n")
	write(paste("time /usr/local/R/3.0.3/bin/Rscript ./01.script/NET.Cluster2Module.R", sp, sep=" "),
	      file = scriptFile, append = TRUE, sep = "\n")
	write(paste("time /usr/local/R/3.0.3/bin/Rscript ./01.script/NET.ComputeEigengene.R", sp, sep=" "),
	      file = scriptFile, append = TRUE, sep = "\n")
	cat("\n..........  Submitted Jobs: ..........\n", sep="")
	cat(paste("qsub -q ", queue, " -e ./01.script/e.clusterScript.", sp, " -o ./01.script/o.clusterScript.", sp, " ", scriptFile, sep=""), sep="\n\n")
	system(paste("chmod u+x", scriptFile, sep=" "))
	system(paste("qsub -q ", queue, " -e ./01.script/e.clusterScript.", sp, " -o ./01.script/o.clusterScript.", sp, " ", scriptFile, sep=""))
}

###########################################################################################
###  Build script(s) and submit job(s) to the computer queue(s).			###
###########################################################################################

buildScript(species[1], queue="rcc-m128-30d")
buildScript(species[2], queue="rcc-m128-30d")

###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list = ls())
gc()
