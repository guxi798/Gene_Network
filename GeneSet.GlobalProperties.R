###########################################################
###	@Author		Xi Gu				###
###	@Time		June 05, 2014			###
###	@Function	Global Properties		###
###########################################################

###########################################################################################
###  This script computes statistical properties for all gene sets, then plot histo-	###
###  -grams for each of the properties. Such properties include but not limited to	###
###  connnectivity (based on SCC), connectivity (based on similarityMatrix), co-regularity, dive-	###
###  -rgence. Wideness will be added in newer version.					###
###  If an extra file containing gene sets of interest, this script will generate	###
###  for each gene set of interest a set of histograms, with a red vertical lines	###
###  indicating where this gene set is in the global background with all gene sets.	###
###											###
###  Input: 1) file in gmt format containing all gene sets.				###
###	    2) file in gmt format containing gene sets of interest. (optional)		###
###  Note, gmt format is the same as defined in GSEA website.				###
###											###
###  Output: 1) background histograms in emf format.					###
###	     2) histograms for each gene set of interest if provided.			###
###	     3) differential expression test in xls format for gene sets of interest	###
###		if provided.								###
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
if(!require("devEMF")){
	install.packages("devEMF", dependencies = TRUE)
	library(devEMF)
}

###########################################################################################
###  Read in arguments, the first argument is gmt file containing all gene sets,	###
###  the second argument is optional which contains gene sets of interest.		###
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
setFile = args[1] # "Gene.Set.Of.Interest.gmt"
backgrFile = args[2] # "01.GeneSet.Kegg.gmt"
species = args[3] 

DataDir = "00.data"
RDataDir = "02.RData"
ResultDir = "03.result"
FigureDir = "05.graph"

###########################################################################################
###  Load RData objects, which should have been saved from upstream analyses.		###
###########################################################################################

cat("\n\n########## Loading RData Objects ##########\n", sep="")
simFile = paste(".", species, RDataDir, "01.datExpr", "01.datExpr.QC.RData", sep="/")
cat("  ", simFile, "\n", sep="")
load(file = simFile)

simFile = paste(".", species, RDataDir, "02.similarityMatrix", "03.similarityMatrix.RData", sep="/")
cat("  ", simFile, "\n", sep="")
load(file = simFile)

###########################################################################################
###  Compute statistical properties. COR is a symmetric matrix whose entries represent 	###
###  gene pair-wise spearman CORrelation. Each row represents one gene with a length    ###
###  same with the total number of genes. TOM is similar to COR but stores topological  ###
###  overlap matrix computed from upstream analyses. Hence, connectivity for each gene  ###
###  is calculated by taking average of the entire row. cv stands for coefficient of	###
###  variance and is a vector the same length of gene number. cv is an index for 	###
###  divergence. Co-regularity is computed for each gene set along the way.		###
###########################################################################################

cat("\n\n########## Computing Global Properties ##########\n", sep="")
COR = cor(datExpr, method="spearman")					## 2 by 2 matrix
con.COR = apply(COR, 1, sum)/dim(datExpr)[2]				## each row is a gene
con.similarityMatrix = apply(similarityMatrix, 1, sum)/dim(datExpr)[2]	## each row is a gene
cv = apply(datExpr, 2, function(x){sd(x)/mean(x)})			
names(cv) = colnames(datExpr)

###########################################################################################
###  Define functions which will be used later. As R requires declaiming a function	###
###  before usage, necessary functions are put here before any real data operations.	###
###########################################################################################

cat("\n\n########## Defining Functions ##########\n", sep="")
GeneSet.Histogram <- function(backgr.con.COR = NULL,	
			      backgr.con.similarityMatrix = NULL,	
			      backgr.within.COR = NULL,	
			      backgr.cv = NULL,		
			      GeneSet.con.COR = NULL,
			      GeneSet.con.similarityMatrix = NULL,
			      GeneSet.within.COR = NULL,
			      GeneSet.cv = NULL,
			      prefix,
				  path){
##### figure 1, distribution of connectivity based on SCC ######
	if(!is.null(backgr.con.COR)){	
		figureFile = paste(path, paste(prefix, "COR.emf", sep="."), sep="/")
		emf(file = figureFile, width = 6, height = 6, family="Helvetica")
	
		h = hist(backgr.con.COR,						## histogram
			breaks=5, 
			xlim = c(-0.1, 0.1),
			main = "Distribution of Connectivity (SCC)",
			xlab = "Connectivity",
			ylab = "Frequency of Gene Sets")
		xfit = seq(min(backgr.con.COR), max(backgr.con.COR), length=40)		## fit the data
		yfit = dnorm(xfit, mean=mean(backgr.con.COR), sd=sd(backgr.con.COR)) 
		yfit = yfit*diff(h$mids[1:2])*length(backgr.con.COR)			
		lines(xfit, yfit, col="blue", lwd=2)					## plot fitted line

		if(!is.null(GeneSet.con.COR)){
			abline(v=GeneSet.con.COR, col="red", lwd=2)			## add line for set of interest
		}
	
		dev.off()
	}
	
##### figure 2, distribution of connectivity base on similarityMatrix #####
	if(!is.null(backgr.con.similarityMatrix)){
		figureFile = paste(path, paste(prefix, "similarityMatrix.emf", sep="."), sep="/")
		emf(file = figureFile, width = 6, height = 6, family="Helvetica")

		h = hist(backgr.con.similarityMatrix,						## histogram
			breaks=5, 
			xlim = c(0, 0.08),
			main = "Distribution of Connectivity (similarityMatrix)",
			xlab = "Connectivity",
			ylab = "Frequency of Gene Sets")
		xfit = seq(0, max(backgr.con.similarityMatrix), length=40)				## fit the data
		yfit = dnorm(xfit, mean=mean(backgr.con.similarityMatrix), sd=sd(backgr.con.similarityMatrix)) 
		yfit = yfit*diff(h$mids[1:2])*length(backgr.con.similarityMatrix) 
		lines(xfit, yfit, col="blue", lwd=2)					## plot fitted line

		if(!is.null(GeneSet.con.similarityMatrix)){
			abline(v=GeneSet.con.similarityMatrix, col="red", lwd=2)			## add line for set of interest
		}
		
		dev.off()
	}

##### figure 3, distribution of CORrelation within gene set #####
	if(!is.null(backgr.within.COR)){
		figureFile = paste(path, paste(prefix, "CORegularity.emf", sep="."), sep="/")
		emf(file = figureFile, width = 6, height = 6, family="Helvetica")

		h = hist(backgr.within.COR,						## histogram
			breaks = 5, 
			xlim = c(-1, 1),
			main = "Distribution of Co-regularity",
			xlab = "Co-regularity",
			ylab = "Frequency of Gene Sets")
		xfit = seq(min(backgr.within.COR), max(backgr.within.COR), length=40)	## fit the data
		yfit = dnorm(xfit, mean=mean(backgr.within.COR), sd=sd(backgr.within.COR)) 
		yfit = yfit*diff(h$mids[1:2])*length(backgr.within.COR) 
		lines(xfit, yfit, col="blue", lwd=2)					## plot fitted line

		if(!is.null(GeneSet.within.COR)){
			abline(v=GeneSet.within.COR, col="red", lwd=2)			## add line for set of interest
		}
	
		dev.off()
	}

###### figure 4, distribution of cv #####
	if(!is.null(backgr.cv)){
		figureFile = paste(path, paste(prefix, "Divergence.emf", sep="."), sep="/")			
		emf(file = figureFile, width = 6, height = 6, family="Helvetica")

		h = hist(backgr.cv,							## histogram
			breaks = 5, 
			xlim = c(0.5, 3.0),
			main = "Distribution of Divergence",
			xlab = "Divergence",
			ylab = "Frequency of Gene Sets")
		
		if(!is.null(GeneSet.cv)){
			abline(v=GeneSet.cv, col="red", lwd=2)				## add line for set of interest
		}
	
		dev.off()
	}
}

GeneSet.Preprocess <- function(data, datExpr){
	set = unlist(strsplit(data, "\t"))[-c(1,2)]
	set.loc = match(set, colnames(datExpr))
	set = set[!is.na(set.loc)]
	set.loc = set.loc[!is.na(set.loc)]
	set.expr = datExpr[, set.loc]
	set.name = unlist(strsplit(data, "\t"))[2]
	set.name = gsub("-", "", set.name, perl=TRUE)
	set.name = gsub("\\s+", "_", set.name, perl=TRUE)
	return(list(set = set,
		    name = set.name,
		    loc = set.loc,
		    expr = set.expr))
}

###########################################################################################
###  Read in gene sets. For each gene set, then compute connectivity (COR/TOM) and 	###
###  divergence by taking average of the CORresponding values of genes within the set.  ###
###  Co-regularity reflects how genes within the set CORrelate with each other. It is   ###
###  the mean of all possible pair-wise gene SCC (except SCC gene versus itself which	###
###  simply 1). Each property is stored in a vector with the same length as the total   ###
###  number of gene set.								###
###########################################################################################

cat("\n\n########## Computing Background Gene Set Properties ##########\n", sep="")
backgrdata = readLines(paste(".", species, ResultDir, "08.GSEA", backgrFile, sep="/"))

backgr.con.COR = NULL
backgr.con.similarityMatrix = NULL
backgr.within.COR = NULL
backgr.cv = NULL
backgr.name = NULL

for(i in 1:length(backgrdata)){
	GeneSet = GeneSet.Preprocess(backgrdata[i], datExpr)
	
	if(length(GeneSet$loc) == 0){							## the set is empty, hence skip it
		cat("\n..........  ", GeneSet$name, " has no genes matched to this species, hence ignored\n\n", sep="")
		next
	}
	cat("\n", sep="")
	
	backgr.con.COR = c(backgr.con.COR, mean(con.COR[GeneSet$loc]))			## connectivity (SCC)
	backgr.con.similarityMatrix = c(backgr.con.similarityMatrix, mean(con.similarityMatrix[GeneSet$loc]))			## connectivity (similarityMatrix)
	backgr.cv = c(backgr.cv, mean(cv[GeneSet$loc]))					## divergence (cv)
	backgr.name = c(backgr.name, GeneSet$name)					## name for all sets
	
	if(length(GeneSet$loc) == 1){
		within.COR = -2								## CORegularity
	}else{
		within.COR = (sum(as.vector(cor(GeneSet$expr, method="spearman")))-length(GeneSet$loc))/(length(GeneSet$loc)^2 - length(GeneSet$loc))
	}
	backgr.within.COR = c(backgr.within.COR, within.COR)
	
	## wideness
}

names(backgr.con.COR) = backgr.name
names(backgr.con.similarityMatrix) = backgr.name
names(backgr.within.COR) = backgr.name
names(backgr.cv) = backgr.name

prefix = "Backgr.GlobalProperties"
GeneSet.Histogram(backgr.con.COR = backgr.con.COR,					## plot the figure
		  backgr.con.similarityMatrix = backgr.con.similarityMatrix,
		  backgr.within.COR = backgr.within.COR,
		  backgr.cv = backgr.cv,
		  prefix = prefix,
		  path = paste(".", species, FigureDir, "06.networkProperties", sep="/"))

###########################################################################################
###  Read in gene sets of interest if provided. For each gene set, then compute		###
###  connectivity (COR/TOM) and divergence by taking average of the corresponding 	###
###  values of genes within the set. Co-regularity reflects how genes within the set    ###
###  CORrelate with each other. It is the mean of all possible pair-wise gene SCC	###
###  (except SCC gene versus itself which simply 1). Each property is stored in a	###
###  vector with the same length as the total number of gene set.			###
###########################################################################################

cat("\n\n########## Computing Properties for Gene Set of Interest ##########\n", sep="")
data = readLines(paste(".", species, DataDir, setFile, sep="/"))

for(i in 1:length(data)){
	GeneSet = GeneSet.Preprocess(data[i], datExpr)

	if(length(GeneSet$loc) == 0){							## the set is empty, hence skip it
		cat("\n..........  ", GeneSet$name, " has no genes matched to both species, hence ignored\n\n", sep="")
		next
	}
	cat("\n", sep="")

	loc = match(GeneSet$name, backgr.name)						## match set of interst to all sets
	GeneSet.con.COR = backgr.con.COR[loc]						## connectivity (SCC)
	GeneSet.con.similarityMatrix = backgr.con.similarityMatrix[loc]						## connectivity (similarityMatrix)
	GeneSet.within.COR = backgr.within.COR[loc]					## divergence (cv)
	GeneSet.cv = backgr.cv[loc]							## CORegularity

	prefix = paste(GeneSet$name, "GlobalProperties", sep=".")		## plot the figure
	GeneSet.Histogram(backgr.con.COR = backgr.con.COR,	
			  backgr.con.similarityMatrix = backgr.con.similarityMatrix,	
			  backgr.within.COR = backgr.within.COR,	
			  backgr.cv = backgr.cv,			
			  GeneSet.con.COR = GeneSet.con.COR,
			  GeneSet.con.similarityMatrix = GeneSet.con.similarityMatrix,
			  GeneSet.within.COR = GeneSet.within.COR,
			  GeneSet.cv = GeneSet.cv,
			  prefix = prefix,
			  path = paste(".", species, FigureDir, "06.networkProperties", sep="/"))
}

###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list = ls())
gc()
