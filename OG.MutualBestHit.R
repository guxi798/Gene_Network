###########################################################
###	@Author		Xi Gu				###
###	@Time		July 8, 2014			###
###	@Function	Find Mutual Best Hit		###
###########################################################

###########################################################################################
###  This script computes eigengene for each module defined in upstream analysis,	###
###  either based on expression data (datExpr) or Topological Overlap Matrix (TOM).	###
###  This part is necessary because computing eigengenes based on TOM is very		###
###  computational intensive. Running the analysis once and storing results can save	###
###  considerable amount of time.							###
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
if(!require("Matrix")){
	install.packages("Matrix", dependencies = TRUE)
	library(Matrix)
}

HomeDir = "/lustre1/escratch1/guxi798_Jun_18/04.Phylo/PhQ.2014Jan24/07.MCL.OG/14.smallData/"
DataDir = "04.all.blastp/"
ResultDir = "14.backblast/"

infile = paste(HomeDir, DataDir, "blast.out", sep="")
outfile = paste(HomeDir, ResultDir, "blast.MBH.tab", sep="")

Ecutoff = 1e-5

data = read.table(infile, header=FALSE, as.is=TRUE)
data = data[data[,11]<Ecutoff, c(1,2,11)]
data = data[data[,1] != data[,2],]

rename <- function(x){
	seg = unlist(strsplit(x, "[|.]"))
	seg = seg[-length(seg)]
	sp = seg[1]
	seg = seg[-1]
	seg = paste(seg, collapse=".")
	c(sp,seg)
}

query = matrix(unlist(lapply(data[,1], rename)), ncol=2, byrow=TRUE)
query.sp = query[,1]
query.gene = query[,2]
hit = matrix(unlist(lapply(data[,2], rename)), ncol=2, byrow=TRUE)
hit.sp = hit[,1]
hit.gene = hit[,2]

combine = paste(query.gene, hit.sp, sep=".")
unik = !duplicated(combine)

mat.name = unique(union(query.gene[unik], hit.gene[unik]))

p = match(query.gene[unik], mat.name)
q = match(hit.gene[unik], mat.name)

mat = sparseMatrix(p, q, x=rep(1, times=length(p)), dims=c(max(p,q), max(p,q)))

for(i in 1:length(p)){
	if(mat[p[i], q[i]]==1 && mat[q[i], p[i]]==1){
		write(paste(mat.name[c(p[i], q[i])], collapse="\t"), outfile, append=TRUE)
	}
}

rm(list=ls())
gc()
