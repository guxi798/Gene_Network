###########################################################
###	@Author		Xi Gu				###
###	@Time		June 06, 2014			###
###	@Function	Categorical Phenotype File	###
###########################################################

###########################################################################################
###  This script generates phenotype file with categorical labels in cls format, as a	###
###  required input for GSEA. In this case, there are two categories: seeds and non-  	###
###  -seeds.										###
###											###
###  Input: 1) sample information file.							###
###											###
###  Output: 1) phenotype file in cls format.						###
###  Here, cls is a categorical phenotype file.						###
###########################################################################################

###########################################################################################
###  Clean and release the memory before start R					###
###########################################################################################

rm(list=ls())
gc()

###########################################################################################
###  Read in arguments. The first arg is genome annotation file and the second arg is	###
###  database name which can be either Kegg or Cyc which should be consistent with the	###
###  annotation file.									###
###########################################################################################

cat("########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
sampleFile = args[1]						## sample file
species = args[2]
mode = args[3]

DataDir = "00.data"
RDataDir = "02.RData"
ResultDir = "03.result"

###########################################################################################
###  Read in data file which summarizes all samples.				 	###
###########################################################################################

cat("\n\n########## Reading in Data File ##########\n", sep="")
file = paste(".", species, DataDir, sampleFile, sep="/")
data = read.csv(file, sep="\t", header=TRUE, as.is=TRUE)
sample.no = length(data$Seed)

###########################################################################################
###  Write data into phenotype file in continuous or categorical cls format.		 	###
###########################################################################################

cat("\n\n########## Writing Phenotype File ##########\n", sep="")
if(mode == "categorical"){
	file = paste(".", species, ResultDir, "08.GSEA", "03.Phenotype.Categorical.cls", sep="/")
	write(paste(sample.no, "2", "1", sep="\t"), file, sep="\n")
	if(data$Seed[1] == "Y"){
		write(paste("#", "Seed", "Non_Seed", sep="\t"), file, sep="\n", append=TRUE)
	}else{
		write(paste("#", "Non_Seed", "Seed", sep="\t"), file, sep="\n", append=TRUE)
	}
	write(paste(data$Seed, collapse="\t"), file, sep="\n", append=TRUE)
}else if(mode == "continuous"){
	cat("\n\n########## Computing Eigengenes ##########\n", sep="")
	
	similarMethod = args[4]
	cat(".......... Based on Similarity Matrix: ", similarMethod, " ..........\n", sep="")
	modFile = paste(".", species, RDataDir, "04.moduleConversion", "04.networkModule.RData", sep="/")
	cat("  ", modFile, "\n", sep="")
	load(modFile)
	
	MEfile = paste(".", species, RDataDir, "04.moduleConversion", paste("04.moduleEigengene", similarMethod, "RData", sep="."), sep="/")
	cat("  ", MEfile, "\n", sep="")
	load(file = MEfile)
	
	for(i in 1:length(MEs)){
		file = paste(".", species, ResultDir, "08.GSEA", paste("03.Phenotype.Eigengene", names(module)[i], "cls", sep="."), sep="/")
		cat(paste(".", species, ResultDir, "08.GSEA", paste("03.Phenotype.Eigengene", names(module)[i], "cls", sep="."), sep="/"), sep="\n")
		
		write("#numeric", file, sep="\n")
		write(paste("#", names(module)[i],"_Eigengene", sep=""), file, sep="\n", append=TRUE)
		write.table(t(MEs[i]), file, quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, sep="\t")
	}
	
}

###########################################################################################
###  Output eigengenes into required format. Note the values of eigengene should have	###
###  been computed in upstream analysis and stored in RData object.			###
###########################################################################################


###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list=ls())
gc()
