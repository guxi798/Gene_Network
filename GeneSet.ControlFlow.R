###########################################################
###	@Author		Xi Gu									###
###	@Time		Aug 01, 2014							###
###	@Function	Top-level Control						###
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
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")

setFile = args[1] # "Gene.Set.Of.Interest.gmt"
backgrFile = args[2] # "01.GeneSet.Kegg.gmt"
OrthologFile = args[3] # "Pvu_Gma_Ortholog.tab"
annotFile = args[4] # "annotation_info_latest.tab"
species = unlist(strsplit(args[5], ","))

###########################################################################################
###  Contruct functions necessary for future call.					###
###########################################################################################

buildScript <- function(scriptFile,
			scriptName,
			parameters,
			queue = "rcc-30d"){
	write("#!/bin/bash", file = scriptFile, append = FALSE, sep="\n\n")
	write(paste("time /usr/local/R/3.0.3/bin/Rscript ./01.script/", scriptName, " ", paste(parameters, collapse=" "), sep=""),
	      file = scriptFile, append = TRUE, sep = "\n")
	cat("\n..........  Submitted Jobs: ..........\n", sep="")
	cat(paste("qsub -q ", queue, " -e ./01.script/e.", scriptName, " -o ./01.script/o.", scriptName, " ", scriptFile, sep=""), sep="\n\n")
	system(paste("chmod u+x", scriptFile, sep=" "))
	system(paste("qsub -q ", queue, " -e ./01.script/e.", parameters[3], ".", scriptName, " -o ./01.script/o.", parameters[3], ".",  scriptName, " ", scriptFile, sep=""))

} 

###########################################################################################
###  Build script(s) and submit job(s) to the computer queue(s).			###
###########################################################################################

##### run GlobalProperties for two species separately #####
buildScript(scriptFile = paste("./01.script/GeneSet.GlobalProperties", species[1], "sh", sep="."),
	    scriptName = "GeneSet.GlobalProperties.R",
	    parameters = c(setFile, backgrFile, species[1]),
	    queue = "rcc-m128-30d")
buildScript(scriptFile = paste("./01.script/GeneSet.GlobalProperties", species[2], "sh", sep="."),
	    scriptName = "GeneSet.GlobalProperties.R",
	    parameters = c(setFile, backgrFile, species[2]),
	    queue = "rcc-m128-30d")

##### run NetworkProperties #####
buildScript(scriptFile = "./01.script/GeneSet.NetworkProperties.sh",
	    scriptName = "GeneSet.NetworkProperties.R",
	    parameters = c(setFile, OrthologFile, annotFile, species[1], species[2]),
	    queue = "rcc-m128-30d")

##### run PlotHeatmap #####
buildScript(scriptFile = "./01.script/GeneSet.PlotHeatmap.sh",
	    scriptName = "GeneSet.PlotHeatmap.R",
	    parameters = c(setFile, OrthologFile, annotFile, species[1], species[2]),
	    queue = "rcc-m128-30d")

