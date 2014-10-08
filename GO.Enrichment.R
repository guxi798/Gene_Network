###########################################################
###	@Author		Xi Gu				###
###	@Time		June 06, 2014			###
###	@Function	Generate Gene Set File		###
###########################################################

###########################################################################################
###  This script generates gene set files in gmt format for GSEA analysis. The file	###
###  contains multiple gene sets, one row per set. The first column is gene set name,  	###
###  the second column is gene set id and the rest are genes assigned to the gene set.	###
###  Gene to set assignment information can be extracted from genome annotation file.	###
###											###
###  Input: 1) genome annotation file.							###
###	    2) database name, currently Kegg and Cyc are available.			###
###											###
###  Output: 1) gene set file in gmt format.						###
###	     2) statistical summary file which records the size for each gene set	###
###  Note, gmt format is defined by GSEA and required in GSEA analysis.			###
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
if(!require("GOSim")){
	install.packages("GOSim", dependencies = TRUE)
	library(GOSim)
}
if(!require("topGO")){
	install.packages("topGO", dependencies = TRUE)
	library(topGO)
}
if(!require("GOSemSim")){
	install.packages("GOSemSim", dependencies = TRUE)
	library(GOSemSim)
}

source("01.script/myfunctions.R")

###########################################################################################
###  Read in arguments. The first arg is genome annotation file and the second arg is	###
###  database name which can be either Kegg or Cyc which should be consistent with the	###
###  annotation file.
###########################################################################################

cat("########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
annotFile = args[1]						## annotation file
species = args[3]

DataDir = "00.data"
RDataDir = "02.RData"
ResultDir = "03.result"
FigureDir = "05.graph"

###########################################################################################
## Construct GO annotation file in required format. 									###
###########################################################################################

GOfile = paste(".", species, DataDir, "GOAnnotation.txt", sep="/")
buildAnnotation(inFile = paste(".", species, DataDir, annotFile, sep="/"), outFile = GOfile)

###########################################################################################
###  Load RData objects, which should have been saved from upstream analyses.			###
###########################################################################################

cat("\n\n########## Loading RData Objects ##########\n", sep="")

exprFile = paste(".", species, RDataDir, "01.datExpr", "01.datExpr.QC.RData", sep="/")
cat("  ", exprFile, "\n", sep="")
load(file = exprFile)

modFile = paste(".", species, RDataDir, "04.moduleConversion", "04.networkModule.RData", sep="/")
load(modFile)
cat("  ", modFile, "\n", sep="")

###########################################################################################
## Run GO enrichment analysis. 															###
###########################################################################################

GOterm = runGO(file = GOfile, module, colnames(datExpr))
save(GOterm, file = paste("", species, RDataDir, "07.GOenrichment", "07.GOenrichment.Unfiltered.RData", sep="/"))

##### write to result file for visualization #####
for(i in 1:length(module)){
	resultFile = paste(".", species, ResultDir, "07.GOenrichment", paste("07.GOenrichment.Unfiltered", names(module)[i], "REVIGO.txt", sep="."), sep="/")
	target = rbind(GOterm[[i]]$MF, GOterm[[i]]$BP)
	target = rbind(target, GOterm[[i]]$CC)
	write(paste("%GO.ID", "p-value", sep="\t"), resultFile, sep="\n")	## output
	write.table(cbind(target$GO.ID, target$pvalue), resultFile, quote=FALSE, row.names=TRUE, col.names=FALSE, append=TRUE, sep="\t")
	
	resultFile = paste(".", species, ResultDir, "07.GOenrichment", paste("07.GOenrichment.Unfiltered", names(module)[i], "Cytoscape.txt", sep="."), sep="/")
	target = rbind(GOterm[[i]]$MF, GOterm[[i]]$BP)
	target = rbind(target, GOterm[[i]]$CC)
	write(paste("GO.ID", "Description", "p-value", "FDR", sep="\t"), resultFile, sep="\n")	## output
	write.table(cbind(target$GO.ID, target$Description, target$pvalue, target$FDR), resultFile, quote=FALSE, row.names=TRUE, col.names=FALSE, append=TRUE, sep="\t")
}

###########################################################################################
### Filter GO terms.																	###
###########################################################################################

GOfiltered = list()
for(i in 1:length(GOterm)){
	MF=GOfilter(data=GOterm[[i]]$MF, AnntProUp=1, LevelUp=100, LevelDown=0, Pvalue=0.05)
	BP=GOfilter(data=GOterm[[i]]$BP, AnntProUp=1, LevelUp=100, LevelDown=0, Pvalue=0.05)
	CC=GOfilter(data=GOterm[[i]]$CC, AnntProUp=1, LevelUp=100, LevelDown=0, Pvalue=0.05)
	
	MFcut=GOSemSim(data=MF,ontology="MF",cut=1)
	BPcut=GOSemSim(data=BP,ontology="BP",cut=1)
	if(dim(CC)[1] > 1){
		CCcut=GOSemSim(data=CC,ontology="CC",cut=1)
	}else{
		CCcut = CC
	}
	
	GOfiltered = c(GOfiltered, list(list(MF=MFcut, BP=BPcut, CC=CCcut)))
}
names(GOfiltered) = names(module)

save(GOfiltered, file = paste("", species, RDataDir, "07.GOenrichment", "07.GOenrichment.Filtered.RData", sep="/"))

##### write to result file for visualization #####
for(i in 1:length(module)){
	resultFile = paste(".", species, ResultDir, "07.GOenrichment", paste("07.GOenrichment.Filtered", names(module)[i], "REVIGO.txt", sep="."), sep="/")
	target = rbind(GOterm[[i]]$MF, GOterm[[i]]$BP)
	target = rbind(target, GOterm[[i]]$CC)
	write(paste("%GO.ID", "p-value", sep="\t"), resultFile, sep="\n")	## output
	write.table(cbind(GO.ID, pvalue), resultFile, quote=FALSE, row.names=TRUE, col.names=FALSE, append=TRUE, sep="\t")
	
	resultFile = paste(".", species, ResultDir, "07.GOenrichment", paste("07.GOenrichment.Filtered", names(module)[i], "Cytoscape.txt", sep="."), sep="/")
	target = rbind(GOterm[[i]]$MF, GOterm[[i]]$BP)
	target = rbind(target, GOterm[[i]]$CC)
	write(paste("%GO.ID", "p-value", sep="\t"), resultFile, sep="\n")	## output
	write.table(cbind(target$GO.ID, target$Description, target$pvalue, target$FDR), resultFile, quote=FALSE, row.names=TRUE, col.names=FALSE, append=TRUE, sep="\t")
}
