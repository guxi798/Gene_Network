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
###  Read in arguments. The first arg is genome annotation file and the second arg is	###
###  database name which can be either Kegg or Cyc which should be consistent with the	###
###  annotation file.
###########################################################################################

cat("########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
annotFile = args[1]						## annotation file
annotDB = tolower(args[2])				## database: Kegg or Cyc
species = args[3]

DataDir = "00.data"
RDataDir = "02.RData"
ResultDir = "03.result"
FigureDir = "05.graph"

###########################################################################################
###  Load RData objects, which should have been saved from upstream analyses.		###
###########################################################################################

cat("\n\n########## Loading RData Objects ##########\n", sep="")
simFile = paste(".", species, RDataDir, "01.datExpr", "01.datExpr.QC.RData", sep="/")	## expression RData object
load(file = simFile)
cat("  ", simFile, "\n", sep="")
genenames = names(datExpr)

cat("\n\n########## Processing Data and Writing to Output File ##########\n", sep="")
if(annotDB == "cyc"){					## Cyc database
	## read in data, extract and reformat information, write to output file
	data = read.csv(paste(".", species, DataDir, annotFile, sep="/"), sep="\t", as.is=TRUE)
	annot = data[, c(1,2,6,8)]				## 1:set id, 2:set name, 6:protein name, 9:gene id
	annot$Gene-name = apply(annot$Gene-name, function(x){	## remove transcript copy id
					temp = unlist(strsplit(x, "[.]")) 
					temp[1]
				})
	gene = unique(annot$Gene-name)				## multiple transcripts corresponding to the same gene
	loc = match(gene, annot$Gene-name)			
	annotUnique = annot[loc,]				## only retain one gene copy
	pathway = unique(annotUnique$Pathway-id)
	
	## write data into output file
	sink(paste(".", species, ResultDir, "08.GSEA", paste("01.GeneSet", annotDB, "gmt", sep="."), sep="/"))
	for(path in pathway){					## write gene set info into file
		isIn = annotUnique$Pathway-id %in% path
		pathName = annotUnique$Pathway-name[isIn][1]	## gene set name
		pathGene = annotUnique$Gene-name[isIn]		## genes assigned to the set
		cat(pathName, path, pathGene, sep="\t")
		cat("\n")
	}
	sink()
	
	## write statistical summary for gene set
	StatFile = paste(".", species, ResultDir, "08.GSEA", paste("05.GeneSet.Module.MatchStatistics", annotDB, "txt", sep="."), sep="/")
	pathName = NULL
	pathID = NULL
	GeneNo = NULL
	for(path in pathway){
		isIn = annotUnique$Pathway-id %in% path
		pathName = c(pathName, annotUnique$Pathway-name[isIn][1]) ## gene set name
		pathGene = annotUnique$Gene-name[isIn]		## genes assigned to the set
		GeneNo = c(GeneNo, sum(pathGene %in% genenames))## number of the genes assigned
		pathID = c(pathID, path)			## gene set id
	}
	pathData = cbind(pathName, pathID, GeneNo)		## combine all info
	pathData = pathData[order(GeneNo, decreasing=TRUE),]	## sort from largest set to smallest set
	write.table(pathData, file=StatFile, row.names=FALSE, quote=FALSE, sep="\t")
								## write into statistical summary file
}else if(annotDB == "kegg"){				## Kegg database
	## read in data, extract and reformat information, write to output file
	data = read.csv(paste(".", species, DataDir, annotFile, sep="/"), sep="\t", as.is=TRUE)	 ##annotation file
	data2 = read.csv(paste(".", species, DataDir, "pathwayname.tab", sep="/"), sep="\t", as.is=TRUE)
								## match annotation info into Kegg pathway
	subData = data[, c(2,9,12)] #2:gene id, 9:kegg id, 12:gene function

	loc = match(subData[,2], data2[,1])
	annot = cbind(subData[!is.na(loc),], data2[loc[!is.na(loc)],])
	## 1:gene, 2:K_id, 3:gene_name, 4:K_id, 5:ko_id, 6:pathname

	gene = unique(annot[,1])
	loc = match(gene, annot[,1])
	annotUnique = annot[loc,]
	pathway = unique(annotUnique[,5])

	## write data into output file
	sink(paste(".", species, ResultDir, "08.GSEA", paste("01.GeneSet", annotDB, "gmt", sep="."), sep="/"))
	for(path in pathway){
		isIn = annotUnique[,5] %in% path
		pathName = annotUnique[,6][isIn][1]		## gene set name
		pathGene = annotUnique[,1][isIn]		## genes assigned to the set
		cat(pathName, path, pathGene, sep="\t")
		cat("\n")
	}
	sink()

	## write statistical summary for gene set
	StatFile = paste(".", species, ResultDir, "08.GSEA", paste("05.GeneSet.Module.MatchStatistics", annotDB, "txt", sep="."), sep="/")
	pathName = NULL
	pathID = NULL
	GeneNo = NULL
	for(path in pathway){
		isIn = annotUnique[,5] %in% path
		pathName = c(pathName, annotUnique[,6][isIn][1])## gene set name
		pathGene = annotUnique[,1][isIn]		## genes assigned to the set
		GeneNo = c(GeneNo, sum(pathGene %in% genenames))## number of the genes assigned
		pathID = c(pathID, path)			## gene set id
	}
	pathData = cbind(pathName, pathID, GeneNo)		## combine all info
	pathData = pathData[order(GeneNo, decreasing=TRUE),]	## sort from largest set to smallest set
	write.table(pathData, file=StatFile, row.names=FALSE, quote=FALSE, sep="\t")
}

##### create soft link point to the gene set gmt file #####
system(paste("ln -sf",  
		paste(".", species, ResultDir, "08.GSEA", paste("01.GeneSet", annotDB, "gmt", sep="."), sep="/"), 
		paste(".", species, DataDir, paste("01.GeneSet", annotDB, "gmt", sep="."), sep="/"),
		sep=" "))

###########################################################################################
###  Clean and release the memory after finish R					###
###########################################################################################

rm(list=ls())
gc()

