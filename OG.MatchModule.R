###########################################################
###	@Author		Xi Gu				###
###	@Time		June 13, 2014			###
###	@Function	Match Modules of Diff Species	###
###########################################################

###########################################################################################
###  This script computes eigengene for each module defined in upstream analysis,	###
###											###
###  Input: 1) expression data stored in RData object.					###
###	    2) gene modules stored in RData object.					###
###	    3) TOM matrix stored in RData object.					###
###											###
###  Output: 1) Module eigengene and ordered eigengenes stored in RData object.		###
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
###  Read in arguments. The first arg is ortholog group file generated from OrthoMCL. 	###
###  The second arg is home directory where data of different species is located. The	###
###  third arg is the head name which will be used for the file name of figures.	###
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
OrthologFile = args[1] #"Pvu_Gma_Ortholog.tab"
species = unlist(strsplit(args[2], ","))

PvuDir = species[1]
GmaDir = species[2]

DataDir = "00.data"
RDataDir = "02.RData"
ResultDir = "03.result"
FigureDir = "05.graph"

###########################################################################################
###  Load RData objects and process ortholog group file, which should have been saved	###
###  from upstream analyses.								###
###########################################################################################

cat("\n\n########## Loading RData Objects ##########\n", sep="")
modFile = paste(".", PvuDir, RDataDir, "04.moduleConversion", "04.networkModule.RData", sep="/")
cat("  ", modFile, "\n", sep="")
load(modFile)

Pvu.module = module

modFile = paste(".", GmaDir, RDataDir, "04.moduleConversion", "04.networkModule.RData", sep="/")
cat("  ", modFile, "\n", sep="")
load(modFile)

Gma.module = module

cat("\n\n########## Reading in Ortholog Groups ##########\n", sep="")
data = readLines(paste(".", PvuDir, DataDir, OrthologFile, sep="/"))
data = data[-1]							## remove head line

###########################################################################################
###  Define necessary functions which will be used in later analysis. Since R is slow	###
###  in for loop, apply pre-defined functions to data object (matrix/vector) should	###
###  increase computation speed.							###
###########################################################################################

cat("\n\n########## Setting Necessary Functions ##########\n", sep="")
Pvu.module.check <- function(x){				## assign OG genes to Phaseolus modules
	group = unlist(strsplit(x, "\t"))
	group = group[grep("Phvul", group)]
	
	judge = NULL
	for(j in 1:length(Pvu.module)){
		judge = c(judge, sum(group %in% Pvu.module[[j]]))
	}
	c(judge, sum(judge), length(group))			## number of genes matched to each module
												## total number of matched genes
												## number of Phaseolus genes in the OG
}

Gma.module.check <- function(x){				## assign OG genes to Glycine modules
	group = unlist(strsplit(x, "\t"))
	group = group[grep("Glyma", group)]
	
	judge = NULL
	for(j in 1:length(Gma.module)){
		judge = c(judge, sum(group %in% Gma.module[[j]]))
	}
	c(judge, sum(judge), length(group))			## number of genes matched to each module
								## total number of matched genes
								## number of Phaseolus genes in the OG
}

Get.Group.Name <- function(x){					## get the name of each OG
	unlist(strsplit(x, "\t"))[1]
}

Compute.Entropy <- function(x){					## compute matching entropy for each OG
	match.gene.no = x[length(x) - 1]
	length(x) = length(x) - 2				## remove last two element
	p = x[x!=0]/match.gene.no
	sum(-p*log(p))/log(match.gene.no)
}

###########################################################################################
###  Apply predefined function to each species.						###
###########################################################################################

cat("\n\n########## Starting Matching Ortholog Groups with Modules ##########\n", sep="")
## apply OG to module assignment function to Phaseolus
Pvu.match = matrix(unlist(lapply(data, Pvu.module.check)), ncol=length(Pvu.module)+2, byrow=TRUE)
colnames(Pvu.match) = c(names(Pvu.module), "match.gene.no", "gene.no")
rownames(Pvu.match) = unlist(lapply(data, Get.Group.Name))

## apply OG to module assignment function to Glycine
Gma.match = matrix(unlist(lapply(data, Gma.module.check)), ncol=length(Gma.module)+2, byrow=TRUE)
colnames(Gma.match) = c(names(Gma.module), "match.gene.no", "gene.no")
rownames(Gma.match) = unlist(lapply(data, Get.Group.Name))

###########################################################################################
###  Summarize some global statistical properties to compare between species. Such	###
###  properties include number of modules each OG is assigned to, and the entropy for	###
###  each module.									###
###########################################################################################

cat("\n\n########## Computing Statistics ##########\n", sep="")
Pvu.match.module = apply(Pvu.match[, 1:length(Pvu.module)]!=0, 1, sum)
Gma.match.module = apply(Gma.match[, 1:length(Gma.module)]!=0, 1, sum)

Pvu.entropy = apply(Pvu.match, 1, Compute.Entropy)
Gma.entropy = apply(Gma.match, 1, Compute.Entropy)

cat("\n\n########## Generating Statistics Reports ##########\n", sep="")
Pvu.report.file = paste(".", PvuDir, ResultDir, "05.moduleMatch", "OrthologGroup2Module.tab", sep="/")
Gma.report.file = paste(".", GmaDir, ResultDir, "05.moduleMatch", "OrthologGroup2Module.tab", sep="/")

cat("The report is ", Pvu.report.file, sep="")
write(paste(paste(names(Pvu.module), collapse="\t"), "module_no", "entropy", sep="\t"), file=Pvu.report.file, append=FALSE, sep="\n")
write.table(cbind(Pvu.match, match.module=Pvu.match.module, entropy=Pvu.entropy), file=Pvu.report.file, append=TRUE)

cat("The report is ", Gma.report.file, sep="")
write(paste(paste(names(Gma.module), collapse="\t"), "module_no", "entropy", sep="\t"), file=Gma.report.file, append=FALSE, sep="\n")
write.table(cbind(Gma.match, match.module=Gma.match.module, entropy=Gma.entropy), file=Gma.report.file)

###########################################################################################
###  Plot distribution of the statistical properties for each species, then plot the 	###
###  comparisons between the species.							###
###########################################################################################

cat("\n\n########## Plotting Statistics Figures ##########\n", sep="")
Pvu.path = paste(".", PvuDir, FigureDir, "05.moduleMatch", sep="/")
Gma.path = paste(".", GmaDir, FigureDir, "05.moduleMatch", sep="/")
cat("The figures are under ", Pvu.path, "****", sep="")

## Distribution.MatchedModule
exp.factor = 1.5
emf(file=paste(Pvu.path, "Distribution.MatchedModule.emf", sep="/"), width=7, height=7)
par(mar = c(5,5,2,2))
barplot(table(Pvu.match.module), 
	cex.lab = exp.factor, 
	cex.axis = exp.factor, 
	cex.names = exp.factor, 
	main = "", 
	xlab = "Number of Matched Module", 
	ylab = "Frequency of Ortholog Groups")
dev.off()

emf(file=paste(Gma.path, "Distribution.MatchedModule.emf", sep="/"), width=7, height=7)
par(mar = c(5,5,2,2))
barplot(table(Gma.match.module), 
	cex.lab = exp.factor, 
	cex.axis = exp.factor, 
	cex.names = exp.factor,
	main = "", 
	xlab = "Number of Matched Module", 
	ylab = "Frequency of Ortholog Groups")
dev.off()

## Distribution.Entropy
exp.factor = 1.5
emf(file=paste(Pvu.path, "Distribution.Entropy.emf", sep="/"), width=7, height=7)
par(mar = c(5,5,2,2))
hist(Pvu.entropy, 
     main = "", 
     cex.lab = exp.factor, 
     cex.axis = exp.factor, 
     xlab = "Entropy", 
     ylab = "Frequency of Ortholog Groups")
dev.off()

emf(file=paste(Gma.path, "Distribution.Entropy.emf", sep="/"), width=7, height=7)
par(mar = c(5,5,2,2))
hist(Gma.entropy, 
     main = "", 
     cex.lab = exp.factor, 
     cex.axis = exp.factor, 
     xlab = "Entropy", 
     ylab = "Frequency of Ortholog Groups")
dev.off()

## Comparison.MatchedModule
set.seed(12345)
noise1 = runif(length(Pvu.match.module), -0.1, 0.1)		## add some noise to the data to distinguish
set.seed(54321)							## overlapped data points
noise2 = runif(length(Pvu.match.module), -0.1, 0.1)

exp.factor = 1.5
emf(file=paste(Pvu.path, "Comparison.MatchedModule.emf", sep="/"), width=7, height=7)
par(mar = c(5,5,2,2))
plot(x = Pvu.match.module + noise1,
     y = Gma.match.module + noise2, 
     type = "p", 
     cex.lab = exp.factor, 
     cex.axis = exp.factor, 
     xlab = "Number of Matched Modules in Phaseolus", 
     ylab = "Number of Matched Modules in Glycine")
model = lm(Gma.match.module~Pvu.match.module)			## add linear regression line
	abline(model, col="red", lwd=exp.factor*2)
	text(7, 4, paste("R.square = ", round(summary(model)$r.squared, digit=4), sep=""), col="red")
dev.off()

## Comparison.Entropy
exp.factor = 1.5
emf(file=paste(Pvu.path, "Comparison.Entropy.emf", sep="/"), width=7, height=7)
par(mar=c(5,5,2,2))
plot(x = Pvu.entropy + noise1, 
     y = Gma.entropy + noise2, 
     type = "p", 
     cex.lab = exp.factor, 
     cex.axis = exp.factor, 
     xlab = "Entropy in Phaseolus", 
     ylab = "Entropy in Glycine")
dev.off()

###########################################################################################
###  This part computes Jaccord Index, Dice's coefficient and Ochiai coefficient.	###
###  Then output results to several files.						###
###########################################################################################

cat("\n\n########## Matching Modules in Different Species ##########\n", sep="")
## reformat OG vs. module relationship for each species ##
Pvu.assign = Pvu.module
Pvu.OG = list()
Pvu.single = NULL
for(i in 1:length(Pvu.assign)){					## do for each module
	names = NULL						## store names of OGs matched to this module
	for(j in 1:length(data)){
		group = unlist(strsplit(data[j], "\t"))
		name = group[1]
		group = group[grep("Phvul", group)]		## get Phaseolus genes
		loc = Pvu.assign[[i]] %in% group
		Pvu.assign[[i]][loc] = paste(name, 1:sum(loc), sep=".")## give matched OG gene an OG name plus an index
		if(sum(loc)>0){names = c(names, name)}
	}
	Pvu.single = c(Pvu.single, length(grep("Phvul", Pvu.assign[[i]])))## number of gene unmapped by any OG
	Pvu.OG = c(Pvu.OG, list(names))				## a list of OG names mapped to this module
}

Gma.assign = Gma.module
Gma.OG = list()
Gma.single = NULL
for(i in 1:length(Gma.assign)){					## do for each module
	names = NULL						## store names of OGs matched to this module
	for(j in 1:length(data)){
		group = unlist(strsplit(data[j], "\t"))
		name = group[1]
		group = group[grep("Glyma", group)]		## get Glycine genes
		loc = Gma.assign[[i]] %in% group
		Gma.assign[[i]][loc] = paste(name, 1:sum(loc), sep=".")## give matched OG gene an OG name plus an index
		if(sum(loc)>0){names = c(names, name)}
	}
	Gma.single = c(Gma.single, length(grep("Glyma", Gma.assign[[i]])))## number of gene unmapped by any OG
	Gma.OG = c(Gma.OG, list(names))				## a list of OG names mapped to this module
}

## compute similarity index at gene level ##
## note since it's unfair to compare species at gene level due to different genome size, this part will not ##
## be reported to an output file, but just offer an extra option/angle to analyze data.			    ##
Jaccord.gene = NULL
Dice.gene = NULL
Ochiai.gene = NULL
for(i in 1:length(Pvu.assign)){
	for(j in 1:length(Gma.assign)){
		Pvu = Pvu.assign[[i]][grep("mcl", Pvu.assign[[i]])]
		Gma = Gma.assign[[j]][grep("mcl", Gma.assign[[j]])]
		n.intersect = length(intersect(Pvu, Gma))
		n.union = length(union(Pvu, Gma))
		n.Pvu = length(Pvu)
		n.Gma = length(Gma)
		Jaccord.gene = c(Jaccord.gene, n.intersect/n.union)		## Jaccord index
		Dice.gene = c(Dice.gene, n.intersect/(n.Pvu + n.Gma))		## Dice's coefficient
		Ochiai.gene = c(Ochiai.gene, n.intersect/sqrt(n.Pvu * n.Gma))	## Ochiai coefficient
	}
}
Jaccord.gene = matrix(Jaccord.gene, ncol=length(Gma.assign), byrow=TRUE)		## change to matrix
	Dice.gene = matrix(Dice.gene, ncol=length(Gma.assign), byrow=TRUE)		
	Ochiai.gene = matrix(Ochiai.gene, ncol=length(Gma.assign), byrow=TRUE)		
colnames(Jaccord.gene) = paste("Gma", names(Gma.module), sep=".")			## rename the matrix
rownames(Jaccord.gene) = paste("Pvu", names(Pvu.module), sep=".")
	colnames(Dice.gene) = paste("Gma", names(Gma.module), sep=".")			
	rownames(Dice.gene) = paste("Pvu", names(Pvu.module), sep=".")
	colnames(Ochiai.gene) = paste("Gma", names(Gma.module), sep=".")			
	rownames(Ochiai.gene) = paste("Pvu", names(Pvu.module), sep=".")
Jaccord.gene = round(Jaccord.gene, digit=4)						## round to 4 digits
	Dice.gene = round(Dice.gene, digit=4)						
	Ochiai.gene = round(Ochiai.gene, digit=4)						

## compute similarity index at ortholog group level ##
Jaccord.OG = NULL
Dice.OG = NULL
Ochiai.OG = NULL
for(i in 1:length(Pvu.OG)){
	for(j in 1:length(Gma.OG)){
		n.intersect = length(intersect(Pvu.OG[[i]], Gma.OG[[j]]))
		n.union = length(union(Pvu.OG[[i]], Gma.OG[[j]]))
		n.Pvu = length(Pvu.OG[[i]])
		n.Gma = length(Gma.OG[[j]])
		Jaccord.OG = c(Jaccord.OG, n.intersect/n.union) 			## Jaccord coefficient
		Dice.OG = c(Dice.OG, n.intersect/(n.Pvu + n.Gma))			## Dice's coefficient
		Ochiai.OG = c(Ochiai.OG, n.intersect/sqrt(n.Pvu * n.Gma))		## Ochiai coefficient
	}
}
Jaccord.OG = matrix(Jaccord.OG, ncol=length(Gma.assign), byrow=TRUE)			## change to matrix
	Dice.OG = matrix(Dice.OG, ncol=length(Gma.assign), byrow=TRUE)		
	Ochiai.OG = matrix(Ochiai.OG, ncol=length(Gma.assign), byrow=TRUE)		
colnames(Jaccord.OG) = paste("Gma", names(Gma.module), sep=".")				## rename the matrix
rownames(Jaccord.OG) = paste("Pvu", names(Pvu.module), sep=".")
	colnames(Dice.OG) = paste("Gma", names(Gma.module), sep=".")			
	rownames(Dice.OG) = paste("Pvu", names(Pvu.module), sep=".")
	colnames(Ochiai.OG) = paste("Gma", names(Gma.module), sep=".")			
	rownames(Ochiai.OG) = paste("Pvu", names(Pvu.module), sep=".")
Jaccord.OG = round(Jaccord.OG, digit=4)							## round to 4 digits
	Dice.OG = round(Dice.OG, digit=4)						
	Ochiai.OG = round(Ochiai.OG, digit=4)						

cat("\n\n########## Writting Similarity Index to Output File ##########\n", sep="")
matchFile = paste(".", PvuDir, ResultDir, "05.moduleMatch", "PvuMod2GmaMod.Jaccord.Gene.xls", sep="/")
	write.table(Jaccord.gene, file=matchFile, quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")
matchFile = paste(".", PvuDir, ResultDir, "05.moduleMatch", "PvuMod2GmaMod.Dice.Gene.xls", sep="/")
	write.table(Dice.gene, file=matchFile, quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")
matchFile = paste(".", PvuDir, ResultDir, "05.moduleMatch", "PvuMod2GmaMod.Ochiai.Gene.xls", sep="/")
	write.table(Ochiai.gene, file=matchFile, quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

matchFile = paste(".", PvuDir, ResultDir, "05.moduleMatch", "PvuMod2GmaMod.Jaccord.OG.xls", sep="/")
	write.table(Jaccord.OG, file=matchFile, quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")
matchFile = paste(".", PvuDir, ResultDir, "05.moduleMatch", "PvuMod2GmaMod.Dice.OG.xls", sep="/")
	write.table(Dice.OG, file=matchFile, quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")
matchFile = paste(".", PvuDir, ResultDir, "05.moduleMatch", "PvuMod2GmaMod.Ochiai.OG.xls", sep="/")
	write.table(Ochiai.OG, file=matchFile, quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

###########################################################################################
###  Second way to match module is to 	###
###########################################################################################

cat("\n\n########## Matching Modules in Different Species ##########\n", sep="")
cat("\n\n.......... Based on Jaccord Coefficients ............\n", sep="")

#loc = (Pvu.match.module<3) & (Pvu.match.module>0) & (Gma.match.module<3) & (Gma.match.module>0)
#match.module = cbind(Pvu.match.module, Gma.match.module)[loc,]
#table(match.module[,1], match.module[,2])

loc = (Pvu.match.module==1) & (Gma.match.module==1)
Pvu.module.max = names(Pvu.module)[apply(Pvu.match[,1:length(Pvu.module)], 1, which.max)[loc]]
Gma.module.max = names(Gma.module)[apply(Gma.match[,1:length(Gma.module)], 1, which.max)[loc]]
Relation = table(Pvu.module.max, Gma.module.max)
cat("\n\n########## Summarizing Results ##########\n", sep="")
table(Relation)

###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list=ls())
gc()
