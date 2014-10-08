###########################################################
###	@Author		Xi Gu				###
###	@Time		Aug 01, 2014			###
###	@Function	Plot Heatmap for Expression	###
###########################################################

###########################################################################################
###  This script computes statistical properties for all gene sets, then plot histo-	###
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

###########################################################################################
###  Read in arguments, the first argument is gmt file containing all gene sets,	###
###  the second argument is optional which contains gene sets of interest.		###
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
setFile = args[1] # "Gene.Set.Of.Interest.gmt"
OrthologFile = args[2] # "Pvu_Gma_Ortholog.tab"
annotFile = args[3] # "annotation_info_latest.tab"
PvuDir = args[4]
GmaDir = args[5]

tissueFile = "record.txt"
DataDir = "00.data"
RDataDir = "02.RData"
ResultDir = "03.result"
FigureDir = "05.graph"

###########################################################################################
###  Load RData objects, which should have been saved from upstream analyses.		###
###########################################################################################

cat("\n\n########## Loading RData Objects ##########\n", sep="")
simFile = paste(".", PvuDir, RDataDir, "01.datExpr", "01.datExpr.QC.RData", sep="/")
cat("  ", simFile, "\n", sep="")
load(file = simFile)

Pvu.datExpr = datExpr[c(1:9, 12, 13, 15:26), ]
Pvu.datExpr = Pvu.datExpr[c(10,13,14,17,15,16,20,18,11,2,21,22,1,12,23,7,6,4,5,3,8,9), ]
Pvu.rank = (dim(Pvu.datExpr)[2] - t(apply(Pvu.datExpr, 1, rank)) + 1)/dim(Pvu.datExpr)[2]
colnames(Pvu.rank) = colnames(Pvu.datExpr)

modFile = paste(".", PvuDir, RDataDir, "04.moduleConversion", "04.networkModule.RData", sep="/")
cat("  ", modFile, "\n", sep="")
load(modFile)
Pvu.modColors = modColors

simFile = paste(".", GmaDir, RDataDir, "04.moduleConversion", "04.networkModule.RData", sep="/")
cat("  ", simFile, "\n", sep="")
load(file = simFile)
Gma.datExpr = datExpr[c(3,7:24), ]
Gma.datExpr = Gma.datExpr[c(2,3,4,5,6,1,13,14,12,17,15,18,16,19,7,11,9,8,10), ]
Gma.rank = (dim(Gma.datExpr)[2] - t(apply(Gma.datExpr, 1, rank)) + 1)/dim(Gma.datExpr)[2]
colnames(Gma.rank) = colnames(Gma.datExpr)

modFile = paste(".", GmaDir, RDataDir, "04.moduleConversion", "04.networkModule.RData", sep="/")
cat("  ", modFile, "\n", sep="")
load(modFile)
Gma.modColors = modColors

###########################################################################################
###  Read ortholog file, gene set file and genome annotation file. The ortholog info	###
###  will be converted into ortholog group object.					###
###########################################################################################

cat("\n\n########## Reading in Gene Set of Interest ##########\n", sep="")
Pvu.tissue = read.csv(paste(".", PvuDir, DataDir, tissueFile, sep="/"), sep="\t")
Pvu.tissue = Pvu.tissue[c(1:9, 12, 13, 15:26), ]
Pvu.tissue = Pvu.tissue[c(10,13,14,17,15,16,20,18,11,2,21,22,1,12,23,7,6,4,5,3,8,9),]
Pvu.data = readLines(paste(".", PvuDir, DataDir, setFile, sep="/"))

Gma.tissue = read.csv(paste(".", GmaDir, DataDir, tissueFile, sep="/"), sep="\t")
Gma.tissue = Gma.tissue[c(3,7:24), ]
Gma.tissue = Gma.tissue[c(2,3,4,5,6,1,13,14,12,17,15,18,16,19,7,11,9,8,10),]
Gma.data = readLines(paste(".", GmaDir, DataDir, setFile, sep="/"))

cat("\n\n########## Reading in Ortholog Groups ##########\n", sep="")
data.OG = readLines(paste(".", PvuDir, DataDir, OrthologFile, sep="/"))
data.OG = data.OG[-1]

groups = list()
groups.name = NULL
for(i in 1:length(data.OG)){
	group = unlist(strsplit(data.OG[i], "\t"))
	groups.name = c(groups.name, group[1])
	groups = c(groups, list(group[-1]))
}
names(groups) = groups.name

cat("\n\n########## Writing Differential Expression Test Files ##########\n", sep="")
Pvu.statFile = paste(".", PvuDir, ResultDir, "06.networkProperties", "Differential.Expression.Test.xls", sep="/")
Gma.statFile = paste(".", GmaDir, ResultDir, "06.networkProperties", "Differential.Expression.Test.xls", sep="/")

write("# This file contains pvalue for differential expression test based on student t test and Mann-Whitney-Wilcoxon Test", file=Pvu.statFile, sep="\n")
write("# This file contains pvalue for differential expression test based on student t test and Mann-Whitney-Wilcoxon Test", file=Gma.statFile, sep="\n")

###########################################################################################
###  Define functions necessary for downstream analysis.				###
###########################################################################################

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

GeneSet.Heatmap <- function(GeneSet, 
			    file, 
			    tissue,
			    modColors,
			    col.extra, 
			    col, 
			    width = 6, 
			    height = 6,
			    colsep = 7){
	emf(file = file, width = width, height = height, family="Helvetica")

	ColSideColors = ifelse(tissue$Seed=="Y", "gray", "palegreen")
	RowSideColors = modColors[GeneSet$loc]

	rowsep = NULL
	for(j in 2:length(col.extra)){
		if(col.extra[j-1] != col.extra[j]){
			rowsep = c(rowsep, j-1)
		}
	}

	heatmap.2(t(GeneSet$expr), 
		scale = "row",
		Rowv = FALSE,
		Colv = FALSE,
		dendrogram = "none",
		main = "", 
		density.info = "none", 
		trace = "none", 
		keysize = 1.0,
		margin = c(16,15), 
		col = col, 
		rowsep = rowsep,
		colsep = colsep,
		sepcolor = "black",
		cexRow = 1.5,
		cexCol = 1.5,
		srtCol = 50,
		RowSideColors = RowSideColors,
		ColSideColors = ColSideColors
	)
	dev.off()
}

GeneSet.LinePlot <- function(file,
			     Pvu.seed,
			     Pvu.nonseed,
			     Gma.seed,
			     Gma.nonseed,
			     Pvu.ratio.loc,
			     Gma.ratio.loc,
			     Pvu.GeneSet,
			     Gma.GeneSet,
			     width = 6,
			     height = 6){
	## Phaseolus
	if(sum(Pvu.ratio.loc)>1){
		Pvu.nonseed.m = apply(Pvu.GeneSet.expr[Pvu.nonseed,Pvu.ratio.loc], 2, mean)
	}else{
		Pvu.nonseed.m = mean(Pvu.GeneSet.expr[Pvu.nonseed,Pvu.ratio.loc])
	}
	Pvu.seed.ratio = Pvu.GeneSet.expr[Pvu.seed,Pvu.ratio.loc]/Pvu.nonseed.m
	Pvu.rank.loc = match(colnames(Pvu.GeneSet.expr)[Pvu.ratio.loc], colnames(Pvu.datExpr))
	Pvu.seed.rank = Pvu.rank[Pvu.seed, Pvu.rank.loc]

	## Glycine
	Gma.ratio.loc = (Gma.col.extra == division[m])
	if(sum(Gma.ratio.loc)>1){
		Gma.nonseed.m = apply(Gma.GeneSet.expr[Gma.nonseed,Gma.ratio.loc], 2, mean)
	}else{
		Gma.nonseed.m = mean(Gma.GeneSet.expr[Gma.nonseed,Gma.ratio.loc])
	}
	Gma.seed.ratio = Gma.GeneSet.expr[Gma.seed,Gma.ratio.loc]/Gma.nonseed.m
	Gma.rank.loc = match(colnames(Gma.GeneSet.expr)[Gma.ratio.loc], colnames(Gma.datExpr))
	Gma.seed.rank = Gma.rank[Gma.seed, Gma.rank.loc]

	emf(file = file, width = width, height = height, family="Helvetica")
	
	ymax = max(cbind(Pvu.seed.ratio, Gma.seed.ratio))
	cat(" Ratio ymax: ", ymax, "\n", sep="")
	
	plot(x = 1:length(Pvu.seed),
	     y = Pvu.seed,
	     ylim = c(0, ymax),
	     ylab = "Seed to Leaf Ratio",
	     xlab = "",
	     main = "",
	     xaxt = "n",
	     type = "n")
	axis(1, at=1:length(Pvu.seed), labels=c("Globular_Seed", "Heart_Seed", "Cotyledon_Seed", "Early_Mature_Seed", "Dry_Seed"))

	if(sum(Pvu.ratio.loc)>1){
		for(k in 1:dim(Pvu.seed.ratio)[2]){
			points(Pvu.seed.ratio[,k], type = "o", col = "black")
		}
	}else{
			points(Pvu.seed.ratio, type = "o", col = "black")
	}
		
	if(sum(Gma.ratio.loc)>1){
		for(k in 1:dim(Gma.seed.ratio)[2]){
			points(Gma.seed.ratio[,k], type = "o", col = "red")
		}
	}else{
			points(Gma.seed.ratio, type = "o", col = "red")
	}
			
	dev.off()
}

###########################################################################################
###  Plot heatmaps for each pathway with ortholog groups separated. Test differential	###
###  expression for each gene between seeds and leaves. Plot seed to leaf ratio of expr	###
###  for each gene.									###
###########################################################################################

##### loop over each pathway #####
for(i in 1:length(Pvu.data)){

##### read in one pathway, check its availability in the expressed data #####
	Pvu.GeneSet = GeneSet.Preprocess(Pvu.data[i], Pvu.datExpr)
	Gma.GeneSet = GeneSet.Preprocess(Gma.data[i], Gma.datExpr)

	#if(is.na(GeneSet.name)){next}
	cat("\n..........  ", Pvu.GeneSet$name, "\t", Gma.GeneSet$name, "  ...........\n", sep="")
	cat("\n..........   ", length(Pvu.GeneSet$loc), "\t", length(Gma.GeneSet$loc), " \n", sep="")
	if(length(Pvu.GeneSet$loc) == 0 | length(Gma.GeneSet$loc) == 0){
		cat("\n..........  ", Pvu.GeneSet$name, " has no genes matched to both species, hence ignored\n\n", sep="")
		next
	}
	cat("\n", sep="")

##### match pathway genes to ortholog groups ####
	m = 0

	## Phaseolus
	Pvu.row.extra = ifelse(Pvu.tissue$Seed=="Y", 10000000, 20000000)
	names(Pvu.row.extra) = rownames(Pvu.datExpr)
	Pvu.col.extra = rep(m, times = length(Pvu.GeneSet.expr[1,]))
	Pvu.group.match = rep(NA, times = length(Pvu.GeneSet.expr[1,]))
	
	## Glycine
	Gma.row.extra = ifelse(Gma.tissue$Seed=="Y", 10000000, 20000000)
	names(Gma.row.extra) = rownames(Gma.datExpr)
	Gma.col.extra = rep(m, times = length(Gma.GeneSet.expr[1,]))
	Gma.group.match = rep(NA, times = length(Gma.GeneSet.expr[1,]))

	for(j in 1:length(groups)){
		group = groups[[j]]
		Pvu.loc = match(Pvu.GeneSet, group)
		Gma.loc = match(Gma.GeneSet, group)
		if((sum(!is.na(Pvu.loc)) > 0) & (sum(!is.na(Gma.loc)) > 0)){
			m = m + 10000000
			Pvu.col.extra[!is.na(Pvu.loc)] = m
			Gma.col.extra[!is.na(Gma.loc)] = m
			#group.match = c(group.match, groups.name[i])
		}
	}

	# Phaseolus
	Pvu.GeneSet.expr = Pvu.GeneSet.expr[, Pvu.col.extra > 0]
	Pvu.GeneSet.loc = Pvu.GeneSet.loc[Pvu.col.extra > 0]
	Pvu.col.extra = Pvu.col.extra[Pvu.col.extra > 0]
	Pvu.GeneSet.expr = Pvu.GeneSet.expr[, order(Pvu.col.extra)]
	Pvu.GeneSet.loc = Pvu.GeneSet.loc[order(Pvu.col.extra)]
	Pvu.col.extra = Pvu.col.extra[order(Pvu.col.extra)]

	# Glycine
	Gma.GeneSet.expr = Gma.GeneSet.expr[, Gma.col.extra > 0]
	Gma.GeneSet.loc = Gma.GeneSet.loc[Gma.col.extra > 0]
	Gma.col.extra = Gma.col.extra[Gma.col.extra > 0]
	Gma.GeneSet.expr = Gma.GeneSet.expr[, order(Gma.col.extra)]
	Gma.GeneSet.loc = Gma.GeneSet.loc[order(Gma.col.extra)]
	Gma.col.extra = Gma.col.extra[order(Gma.col.extra)]

##### plot heatmap for expression across tissues #####

	cat("\n\n########## Generating Figures for Gene Sets ##########\n", sep="")
	#my.palette = colorRampPalette(c("darkturquoise", "ivory", "orchid1"))(n=299)
	my.palette = colorRampPalette(c("blue", "ivory", "red"))(n=299)
	#layout(matrix(c(1,2), 2, 1, byrow=TRUE))

	##### Phaseolus
	figureFile = paste(".", PvuDir, FigureDir, "06.networkProperties", paste(Pvu.GeneSet.name, "Heatmap.datExpr.emf", sep="."), sep="/")
	GeneSet.Heatmap(GeneSet = Pvu.GeneSet,
			file = figureFile,
			tissue = Pvu.tissue,
			modColors = Pvu.modColors,
			col.extra = Pvu.col.extra,
			col = my.palette,
			width = 10,
			height = 6, 
			colsep = 7)


	#### Glycine ####
	figureFile = paste(".", GmaDir, FigureDir, "06.networkProperties", paste(Pvu.GeneSet.name, "Heatmap.datExpr.emf", sep="."), sep="/")
	GeneSet.Heatmap(GeneSet = Gma.GeneSet,
			file = figureFile,
			tissue = Gma.tissue,
			modColors = Gma.modColors,
			col.extra = Gma.col.extra,
			col = my.palette,
			width = 10,
			height = 9, 
			colsep = 5)
	
##### calculate significance of differential expression ####
	## Phaseolus
	Pvu.pvalue = NULL
	#Pvu.seed = 1:7
	#Pvu.nonseed = 9:15
	Pvu.seed = c(1,2,5:7)
	Pvu.nonseed = c(9,14,15)
	
	for(j in 1:dim(Pvu.GeneSet$expr)[2]){
		expr = Pvu.GeneSet$expr[, j]
		t.two = t.test(expr[Pvu.seed], expr[Pvu.nonseed], alternative="two.side")$p.value
		t.less = t.test(expr[Pvu.seed], expr[Pvu.nonseed], alternative="less")$p.value
		t.greater = t.test(expr[Pvu.seed], expr[Pvu.nonseed], alternative="greater")$p.value
		w.two = wilcox.test(expr[Pvu.seed], expr[Pvu.nonseed], alternative="two.side")$p.value
		w.less = wilcox.test(expr[Pvu.seed], expr[Pvu.nonseed], alternative="less")$p.value
		w.greater = wilcox.test(expr[Pvu.seed], expr[Pvu.nonseed], alternative="greater")$p.value
		Pvu.pvalue = rbind(Pvu.pvalue, c(t.two,t.less,t.greater,w.two,w.less,w.greater))
	}
	rownames(Pvu.pvalue) = colnames(Pvu.GeneSet$expr)
	colnames(Pvu.pvalue) = c("Two Side T", "Less T", "Greater T", "Two Side W", "Less W", "Greater W")
	write(paste("# ", Pvu.GeneSet$name, sep=""), file=Pvu.statFile, append=TRUE, sep="\n")
	write.table(round(Pvu.pvalue, digits=4), file=Pvu.statFile, row.names=TRUE, col.names=TRUE, quote=FALSE, append=TRUE, sep="\t")
	
	## Glycine
	Gma.pvalue = NULL
	#Gma.seed = 1:5
	#Gma.nonseed = c(7,9,11,13,15)
	Gma.seed = c(1:5)
	Gma.nonseed = c(7,9,15)
	for(j in 1:dim(Gma.GeneSet$expr)[2]){
		expr = Gma.GeneSet$expr[, j]
		t.two = t.test(expr[Gma.seed], expr[Gma.nonseed], alternative="two.side")$p.value
		t.less = t.test(expr[Gma.seed], expr[Gma.nonseed], alternative="less")$p.value
		t.greater = t.test(expr[Gma.seed], expr[Gma.nonseed], alternative="greater")$p.value
		w.two = wilcox.test(expr[Gma.seed], expr[Gma.nonseed], alternative="two.side")$p.value
		w.less = wilcox.test(expr[Gma.seed], expr[Gma.nonseed], alternative="less")$p.value
		w.greater = wilcox.test(expr[Gma.seed], expr[Gma.nonseed], alternative="greater")$p.value
		Gma.pvalue = rbind(Gma.pvalue, c(t.two,t.less,t.greater,w.two,w.less,w.greater))
	}
	rownames(Gma.pvalue) = colnames(Gma.GeneSet$expr)
	colnames(Gma.pvalue) = c("Two Side T", "Less T", "Greater T", "Two Side W", "Less W", "Greater W")
	write(paste("# ", Gma.GeneSet$name, sep=""), file=Gma.statFile, append=TRUE, sep="\n")
	write.table(Gma.pvalue, file=Gma.statFile, row.names=TRUE, col.names=TRUE, quote=FALSE, append=TRUE, sep="\t")

##### plot seed to green ratio	#####
	Pvu.seed = c(1,2,5:7)
	Pvu.nonseed = c(9,14,15)
	Gma.seed = c(1:5)
	Gma.nonseed = c(7,9,15)
	
	division = unique(Pvu.col.extra)
	
	for(m in 1:length(division)){
		Pvu.ratio.loc = (Pvu.col.extra == division[m])
		Gma.ratio.loc = (Gma.col.extra == division[m])
					
		figureFile = paste(".", PvuDir, FigureDir, "06.networkProperties", paste(Pvu.GeneSet.name, "Seed2LeafRatio.OG", m, "emf", sep="."), sep="/")
		GeneSet.LinePlot(file = figureFile,
				 Pvu.seed = Pvu.seed,
				 Pvu.nonseed = Pvu.nonseed,
				 Gma.seed = Gma.seed,
				 Gma.nonseed = Gma.nonseed,
				 Pvu.ratio.loc = Pvu.ratio.loc ,
				 Gma.ratio.loc = Gma.ratio.loc,
				 width = 10,
				 height = 6)

	## plot expression rank
		figureFile = paste(".", PvuDir, FigureDir, "06.networkProperties", paste(Pvu.GeneSet.name, ".SeedRank.OG", m, "emf", sep="."), sep="/")
		emf(file = figureFile, 
			width = 10, 
			height = 6, 
			family="Helvetica"
		)
		ymax = max(cbind(Pvu.seed.rank, Gma.seed.rank))
		cat(" Rank ymax: ", ymax, "\n", sep="")
	
		plot(x = 1:length(Pvu.seed),
			y = Pvu.seed,
			ylim = c(0, ymax),
			ylab = "Gene Rank of Expression",
			xlab = "",
			main = "",
			xaxt = "n",
			type = "n")
		axis(1, at=1:length(Pvu.seed), 
			labels=c("Globular_Seed", "Heart_Seed", "Cotyledon_Seed", "Early_Mature_Seed", "Dry_Seed"))

		if(sum(Pvu.ratio.loc)>1){
			for(k in 1:dim(Pvu.seed.rank)[2]){
				points(Pvu.seed.rank[,k], type = "o", col = "black")
			}
		}else{
				points(Pvu.seed.rank, type = "o", col = "black")
		}
		
		if(sum(Gma.ratio.loc)>1){
			for(k in 1:dim(Gma.seed.rank)[2]){
				points(Gma.seed.rank[,k], type = "o", col = "red")
			}
		}else{
				points(Gma.seed.rank, type = "o", col = "red")
		}
			
	dev.off()

	}

}

###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list = ls())
gc()
