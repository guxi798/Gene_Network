###########################################################
###	@Author		Xi Gu				###
###	@Time		June 11, 2014			###
###	@Function	Rank the Gene List		###
###########################################################

###########################################################################################
###  This script generates pre-ordered gene list. GSEA only provides several ways to	###
###  rank the genes, e.g. Pearson Correlation Coefficient, but it also allows users to	###
###  input their own ordered gene list. This script rank the genes based on Spearman	###
###  Correlation Coefficient which is not available in GSEA internal code but is quite	###
###  popular due to its robustness to outliers.						###
###											###
###  Input: 1) expression data stored in RData object.					###
###	    2) gene to module assignment stored in RData object.			###
###	    3) similarity matrix, either "datExpr" or "TOM".				###
###	    4) correlation method, here we use "spearman"				###
###											###
###  Output: 1) ordered gene list file in rnk format.					###
###  Here, rnk is a gene rank file defined by GSEA.					###
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
###  Read in arguments. The first arg is the similarity matrix and the second arg is	###
###  method used to rank gene list. Here only "pearson", "kendall" and "spearman" are	###
###  acceptable.									###
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
#noCluster = args[1]
similarMethod = args[1] # datExpr or similarityMatrix
corMethod = args[2]	# pearson, kendall, spearman
species = args[3]

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
	cat(".......... Loading ", similarMethod, " Object..........\n", sep="")
	load(file = simFile)
	similarityMatrix = datExpr
}else if(similarMethod == "similarityMatrix"){
	simFile = paste(".", species, RDataDir, "02.similarityMatrix", "03.similarityMatrix.RData", sep="/")
	cat(".......... Loading ", similarMethod, " Object..........\n", sep="")
	load(file = simFile)
	similarityMatrix = similarityMatrix
}
cat("  ", simFile, "\n", sep="")

cat(".......... Based on Similarity Matrix: ", similarMethod, " ..........\n", sep="")
MEfile = paste(".", species, RDataDir, "04.moduleConversion", paste("04.moduleEigengene", similarMethod, "RData", sep="."), sep="/")
load(file = MEfile)

###########################################################################################
###  Rank all genes based on the correlation matrix from the input. Output the ordered	###
###  gene list into rnk file (one per module).						###
###########################################################################################

cat("\n\n########## Ranking Genes ##########\n", sep="")
genenames = colnames(similarityMatrix)
for(i in 1:length(module)){
	loc = match(module[[i]], genenames)				## assign genes to modules
	subSimilarityMatrix = similarityMatrix[,loc]
	corList = apply(subSimilarityMatrix, 2, 
			function(x){cor(x, MEs[i], method=corMethod)})	## compute correlation with eigengene
	names(corList) = module[[i]]
	corList = sort(corList, decreasing=TRUE)			## rank genes
	file = paste(".", species, ResultDir, "08.GSEA", paste("04.RankedGeneList.Eigengene", names(module)[i], "rnk", sep="."), sep="/")
	write(paste("# Rank based on Eigengene and ", similarMethod, sep=""), file, sep="\n")	## output
	write.table(as.matrix(corList), file, quote=FALSE, row.names=TRUE, col.names=FALSE, append=TRUE, sep="\t")
}

###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list=ls())
gc()








