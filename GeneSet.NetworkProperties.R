###########################################################
###	@Author		Xi Gu				###
###	@Time		July 29, 2014			###
###	@Function	Network Properties		###
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
if(!require("devEMF")){
	install.packages("devEMF", dependencies = TRUE)
	library(devEMF)
}
if(!require("igraph")){
	install.packages("igraph", dependencies = TRUE)
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
##### Phaseolus #####
dataFile = paste(".", PvuDir, RDataDir, "01.datExpr", "01.datExpr.QC.RData", sep="/")
cat("  ", dataFile, "\n", sep="")
load(file = dataFile)
Pvu.datExpr = datExpr

simFile = paste(".", PvuDir, RDataDir, "02.similarityMatrix", "03.similarityMatrix.RData", sep="/")
cat("  ", simFile, "\n", sep="")
load(file = simFile)
Pvu.similarityMatrix = similarityMatrix

modFile = paste(".", PvuDir, RDataDir, "04.moduleConversion", "04.networkModule.RData", sep="/")
cat("  ", modFile, "\n", sep="")
load(modFile)
Pvu.modColors = modColors

##### Glycine #####
dataFile = paste(".", GmaDir, RDataDir, "01.datExpr", "01.datExpr.QC.RData", sep="/")
cat("  ", dataFile, "\n", sep="")
load(file = dataFile)
Gma.datExpr = datExpr

simFile = paste(".", GmaDir, RDataDir, "02.similarityMatrix", "03.similarityMatrix.RData", sep="/")
cat("  ", simFile, "\n", sep="")
load(file = simFile)
Gma.similarityMatrix = similarityMatrix

modFile = paste(".", GmaDir, RDataDir, "04.moduleConversion", "04.networkModule.RData", sep="/")
cat("  ", modFile, "\n", sep="")
load(modFile)
Gma.modColors = modColors

###########################################################################################
###  Read ortholog file, gene set file and genome annotation file. The ortholog info	###
###  will be converted into ortholog group object.					###
###########################################################################################

Pvu.genenames = colnames(Pvu.datExpr)
Gma.genenames = colnames(Gma.datExpr)

cat("\n\n########## Reading in Gene Set of Interest ##########\n", sep="")
Pvu.data = readLines(paste(".", PvuDir, DataDir, setFile, sep="/"))
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

cat("\n\n########## Reading in Gene Annotation File ##########\n", sep="")
Pvu.annot = read.csv(paste(".", PvuDir, DataDir, annotFile, sep="/"), sep="\t", as.is=TRUE, header=FALSE, na.strings="")
Gma.annot = read.csv(paste(".", GmaDir, DataDir, annotFile, sep="/"), sep="\t", as.is=TRUE, header=FALSE, na.strings="")

Pvu.statFile = paste(".", PvuDir, ResultDir, "06.networkProperties", "GeneSet.NetworkProperties.Summary.txt", sep="/")
Gma.statFile = paste(".", GmaDir, ResultDir, "06.networkProperties", "GeneSet.NetworkProperties.Summary.txt", sep="/")

write("This file records network properties for genes of interest organized by ortholog groups.", file=Pvu.statFile, append=FALSE, sep="\n")
write("This file records network properties for genes of interest organized by ortholog groups.", file=Gma.statFile, append=FALSE, sep="\n")

###########################################################################################
###  Reduce graph size by removing a quantile of weakest edges from similarity matrix.	###
###########################################################################################

##### remove the 25% weakest edges #####
cat("\n\n########## Reducing graph size by removing a quantile of weakest edges ##########\n", sep="")
Pvu.edgeCut = quantile(as.vector(Pvu.similarityMatrix), 0.25)
Pvu.similarityMatrix[Pvu.similarityMatrix <= Pvu.edgeCut] = 0

Gma.edgeCut = quantile(as.vector(Gma.similarityMatrix), 0.25)
Gma.similarityMatrix[Gma.similarityMatrix <= Gma.edgeCut] = 0

###########################################################################################
###  Build graph for each module and store module graphs into a list object for each	###
###  species. This helps reducing graph size as well.					###
###########################################################################################

cat("\n\n########## Building graph for each module ##########\n", sep="")

##### Phaseolus #####
Pvu.graph = list()
Pvu.modNames = unique(Pvu.modColors)
for(Pvu.modName in Pvu.modNames){
	Pvu.loc = which(Pvu.modColors == Pvu.modName)
	Pvu.sim = Pvu.similarityMatrix[Pvu.loc, Pvu.loc]
	Pvu.modGraph = graph.adjacency(Pvu.sim, mode="upper", weight=TRUE, diag=FALSE)
	Pvu.graph = c(Pvu.graph, list(Pvu.modGraph))
}
names(Pvu.graph) = Pvu.modNames
graphFile = paste(".", PvuDir, RDataDir, "06.networkProperties", "06.moduleGraph.RData", sep="/")
save(Pvu.graph, file = graphFile)

##### Glycine #####
Gma.graph = list()
Gma.modNames = unique(Gma.modColors)
for(Gma.modName in Gma.modNames){
	Gma.loc = which(Gma.modColors == Gma.modName)
	Gma.sim = Gma.similarityMatrix[Gma.loc, Gma.loc]
	Gma.modGraph = graph.adjacency(Gma.sim, mode="upper", weight=TRUE, diag=FALSE)
	Gma.graph = c(Gma.graph, list(Gma.modGraph))
}
names(Gma.graph) = Gma.modNames
graphFile = paste(".", GmaDir, RDataDir, "06.networkProperties", "06.moduleGraph.RData", sep="/")
save(Gma.graph, file = graphFile)

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

###########################################################################################
###  Relate ortholog groups between the two species. Calculate network properties,	###
###  including connectivity, betweenness, closeness, average nearest neighbor centrality###
###  and clustering coefficient. Top 10 nearest neighbor genes will also be recorded.	###
###########################################################################################

cat("\n\n########## Relate ortholog groups between the two species ##########\n", sep="")

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

	##### write to the stat file #####
	cat("\n..........  write to the stat file  ...........\n", sep="")
		write(Pvu.GeneSet$name, file=Pvu.statFile, append=TRUE, sep="\n")
	write(Gma.GeneSet$name, file=Gma.statFile, append=TRUE, sep="\n")

##### match pathway genes to ortholog groups #####
	cat("\n..........  match pathway genes to ortholog groups  ...........\n", sep="")
	m = 0

	## Phaseolus
	Pvu.col.extra = rep(m, times = length(Pvu.GeneSet$expr[1,]))
	Pvu.group.match = rep(NA, times = length(Pvu.GeneSet$expr[1,]))
	
	## Glycine
	Gma.col.extra = rep(m, times = length(Gma.GeneSet$expr[1,]))
	Gma.group.match = rep(NA, times = length(Gma.GeneSet$expr[1,]))

	for(j in 1:length(groups)){
		group = groups[[j]]
		Pvu.loc = match(Pvu.GeneSet$set, group)
		Gma.loc = match(Gma.GeneSet$set, group)
		if((sum(!is.na(Pvu.loc)) > 0) & (sum(!is.na(Gma.loc)) > 0)){
			m = m + 1
			Pvu.col.extra[!is.na(Pvu.loc)] = m
			Gma.col.extra[!is.na(Gma.loc)] = m
		}
	}

	# Phaseolus
	Pvu.GeneSet$loc = Pvu.GeneSet$loc[Pvu.col.extra > 0]
	Pvu.col.extra = Pvu.col.extra[Pvu.col.extra > 0]
	Pvu.GeneSet$loc = Pvu.GeneSet$loc[order(Pvu.col.extra)]
	Pvu.col.extra = Pvu.col.extra[order(Pvu.col.extra)]

	# Glycine
	Gma.GeneSet$loc = Gma.GeneSet$loc[Gma.col.extra > 0]
	Gma.col.extra = Gma.col.extra[Gma.col.extra > 0]
	Gma.GeneSet$loc = Gma.GeneSet$loc[order(Gma.col.extra)]
	Gma.col.extra = Gma.col.extra[order(Gma.col.extra)]

##### write file for each pathway #####
	cat("\n..........  write file for each pathway  ...........\n", sep="")
	write(paste("OG_ID", 
		    "Gene_ID", 
		    "Expr_Mean",
		    "Expr_CV",
		    "Expr_SamCoverage",
		    "Connectivity", 
		    #"Betweenness", 
		    "Closeness", 
		    "Average_Nearest_Neigbor_Degree", 
		    "Clustering_Coefficient", 
		    "Neighbor_Genes", 
		    sep="\t"), file=Pvu.statFile, append=TRUE, sep="\n")
	write(paste("OG_ID", 
		    "Gene_ID", 
		    "Expr_Mean",
		    "Expr_CV",
		    "Expr_SamCoverage",
		    "Connectivity", 
		    #"Betweenness",
		    "Closeness",
		    "Average_Nearest_Neigbor_Degree", 
		    "Clustering_Coefficient", 
		    "Neighbor_Genes", 
		    sep="\t"), file=Gma.statFile, append=TRUE, sep="\n")

##### extract network properties and write to file #####
	cat("\n..........  extract network properties and write to file  ...........\n", sep="")
	for(j in 1:m){
		## Phaseolus
		Pvu.gene.locs = Pvu.GeneSet$loc[Pvu.col.extra == j]
		for(Pvu.gene.loc in Pvu.gene.locs){
			Pvu.group = paste("OG", j, sep=".")
			Pvu.gene = Pvu.genenames[Pvu.gene.loc]
			Pvu.expr = Pvu.datExpr[, Pvu.gene.loc]
			Pvu.modLoc = match(Pvu.modColors[Pvu.gene.loc], Pvu.modNames)
			Pvu.G = Pvu.graph[[Pvu.modLoc]]
			Pvu.loc = which(Pvu.modColors == Pvu.modColors[Pvu.gene.loc])
			Pvu.sim = Pvu.similarityMatrix[Pvu.loc, Pvu.loc]

			Pvu.gene.mean = mean(Pvu.expr)							## mean expression
			Pvu.gene.cv = sd(Pvu.expr)/mean(Pvu.expr)					## coefficient variance
			Pvu.gene.sam = sum(Pvu.expr > quantile(as.vector(Pvu.datExpr), 0.05))		## tissue coverage
			Pvu.gene.con = sum(Pvu.sim[match(Pvu.gene, colnames(Pvu.sim)),])		## connectivity
			#Pvu.gene.btn = betweenness(Pvu.G, v=Pvu.gene, directed=FALSE, normalized=TRUE)	## betweenness
			Pvu.gene.cln = closeness(Pvu.G, vids=Pvu.gene, mode="all", normalized=TRUE)	## closeness
			Pvu.gene.knn = graph.knn(Pvu.G, vids=Pvu.gene)					## average nearest neighbor centrality
			Pvu.gene.tra = transitivity(Pvu.G, type="local", vids=Pvu.gene)			## clustering coefficient

			Pvu.neigbor = Pvu.genenames[order(Pvu.sim[match(Pvu.gene, colnames(Pvu.sim)),])[1:11]]
			Pvu.neigbor.annot = Pvu.annot$V13[match(Pvu.neigbor, Pvu.annot$V2)]
			Pvu.neigbor = paste(Pvu.neigbor, Pvu.neigbor.annot, sep=":")
			
			Pvu.info = paste(Pvu.group, 
					 Pvu.gene, 
					 Pvu.gene.mean,
					 Pvu.gene.cv,
					 Pvu.gene.sam,
					 Pvu.gene.con, 
					 #Pvu.gene.btn, 
					 Pvu.gene.cln, 
					 Pvu.gene.knn, 
					 Pvu.gene.tra, 
					 Pvu.neigbor, 
					 sep="\t")
			write(Pvu.info, file= Pvu.statFile, append=TRUE, sep="\t")
		}

		## Glycine
		Gma.gene.locs = Gma.GeneSet$loc[Gma.col.extra == j]
		for(Gma.gene.loc in Gma.gene.locs){
			Gma.group = paste("OG", j, sep=".")
			Gma.gene = Gma.genenames[Gma.gene.loc]
			Gma.modLoc = match(Gma.modColors[Gma.gene.loc], Gma.modNames)
			Gma.G = Gma.graph[[Gma.modLoc]]
			Gma.loc = which(Gma.modColors == Gma.modColors[Gma.gene.loc])
			Gma.sim = Gma.similarityMatrix[Gma.loc, Gma.loc]
			
			Gma.gene.mean = mean(Gma.expr)							## mean expression
			Gma.gene.cv = sd(Gma.expr)/mean(Gma.expr)					## coefficient variance
			Gma.gene.sam = sum(Gma.expr > quantile(as.vector(Gma.datExpr), 0.05))		## tissue coverage
			Gma.gene.con = sum(Gma.sim[match(Gma.gene, colnames(Gma.sim)),])		## connectivity
			#Gma.gene.btn = betweenness(Gma.G, v=Gma.gene, directed=FALSE, normalized=TRUE)	## betweenness
			Gma.gene.cln = closeness(Gma.G, vids=Gma.gene, mode="all", normalized=TRUE)	## closeness
			Gma.gene.knn = graph.knn(Gma.G, vids=Gma.gene)					## average nearest neighbor centrality
			Gma.gene.tra = transitivity(Gma.G, type="local", vids=Gma.gene)			## clustering coefficient

			Gma.neigbor = Gma.genenames[order(Gma.sim[match(Gma.gene, colnames(Gma.sim)),])[1:11]]
			Gma.neigbor.annot = Gma.annot$V13[match(Gma.neigbor, Gma.annot$V2)]
			Gma.neigbor = paste(Gma.neigbor, Gma.neigbor.annot, sep=":")
			
			Gma.info = paste(Gma.group, 
					 Gma.gene, 
					 Gma.gene.mean,
					 Gma.gene.cv,
					 Gma.gene.sam,
					 Gma.gene.con, 
					 #Gma.gene.btn, 
					 Gma.gene.cln, 
					 Gma.gene.knn, 
					 Gma.gene.tra, 
					 Gma.neigbor, 
					 sep="\t")
			write(Gma.info, file= Gma.statFile, append=TRUE, sep="\t")
		}
	}
}


## genes mapped to the network

## genes mapped to the modules


## globally assortative

## locally dissortative

## clustering coefficient



###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list = ls())
gc()


