rm(list=ls())
gc()

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
if(!require("WGCNA")){
	install.packages("WGCNA", dependencies = TRUE)
	library(WGCNA)
}
if(!require("psych")){
	install.packages("psych", dependencies = TRUE)
	library(psych)
}
if(!require("fpc")){
	install.packages("fpc", dependencies = TRUE)
	library(fpc)
}
if(!require("vegan")){
	install.packages("vegan", dependencies = TRUE)
	library(vegan)
}
if(!require("mclust")){
	install.packages("mclust", dependencies = TRUE)
	library(mclust)
}
if(!require("cluster")){
	install.packages("cluster", dependencies = TRUE)
	library(cluster)
}
if(!require("NbClust")){
	install.packages("NbClust", dependencies = TRUE)
	library(NbClust)
}
if(!require("pvclust")){
	install.packages("pvclust", dependencies = TRUE)
	library(pvclust)
}
if(!require("MASS")){
	install.packages("MASS", dependencies = TRUE)
	library(MASS)
}

cat("\n\n########## Reading Arguments ##########\n", sep="")
#args = commandArgs(TRUE)
#for (arg in args) cat("  ", arg, "\n", sep="")

DataDir = "../00.data/"
RDataDir = "../02.RData/01.raw/"
ResultDir = "../03.result/01.raw/03.Match/01.GeneSet/"
FigureDir = "../05.graph/01.raw/04.Cluster/"

cat("\n\n########## Loading RData Objects ##########\n", sep="")
simFile = paste(RDataDir, "01.datExpr.QC.RData", sep="")
	cat("  ", simFile, "\n", sep="")
	load(file = simFile)

simFile = paste(RDataDir, "03.TOM.Unsigned.Signed.RData", sep="")
	cat(".......... Loading TOM Object..........\n", sep="")
	load(file = simFile)

cat("\n\n########## Trying Different Number of Clusters ##########\n", sep="")
## 1. Look for a bend or elbow in the sum of squared error (SSE) screen plot
cat(".......... No.1 elbow in the sum of squared error (SSE) plot ..........\n", sep="")
#figureFile = paste(FigureDir, "03.Kmeans.TOM.ClusterNo.SSE.emf", sep="")
#emf(file = figureFile, 
#	width = 6, 
#	height = 5, 
#	family="Helvetica"
#)

#wss = (nrow(TOM)-1)*sum(apply(TOM, 2, var))
#for(i in 2:15) {
#	wss[i] = sum(kmeans(TOM, centers=i)$withinss)
#}
#plot(1:15, wss, type="b", xlab="Number of Clusters",
#     ylab="Within groups sum of squares")

#dev.off()

## 2. Partition around medoids to estimate the number of clusters
cat(".......... No.2 partition around medoids for estimation ..........\n", sep="")
figureFile = paste(FigureDir, "03.Kmeans.TOM.ClusterNo.PAM.emf", sep="")
emf(file = figureFile, 
	width = 6, 
	height = 10, 
	family="Helvetica"
)

pamk.best = pamk(TOM)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
plot(pam(TOM, pamk.best$nc))

dev.off()

asw = numeric(15)
for (k in 2:15)
  asw[[k]] = pam(TOM, k)$silinfo$avg.width
k.best = which.max(asw)
cat("silhouette-optimal number of clusters:", k.best, "\n")

## 3. Calinsky criterion
cat(".......... No.3 Calinsky criterion ..........\n", sep="")
figureFile = paste(FigureDir, "03.Kmeans.TOM.ClusterNo.Calinsky.emf", sep="")
emf(file = figureFile, 
	width = 6, 
	height = 5, 
	family="Helvetica"
)

fit = cascadeKM(scale(TOM,center = TRUE,scale = TRUE), 1, 15, iter=1000)
plot(fit, 
	sortg = TRUE, 
	grpmts.plot = TRUE)
calinski.best = as.numeric(which.max(fit$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")

## 4. BIC for EM
cat(".......... No.4 Bayesian Information Criterion for Expectation-maximization ..........\n", sep="")
figureFile = paste(FigureDir, "03.Kmeans.TOM.ClusterNo.Calinsky.emf", sep="")
emf(file = figureFile, 
	width = 6, 
	height = 15, 
	family="Helvetica"
)

d_clust = Mclust(as.matrix(TOM), G=1:15)
m.best = dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
plot(d_clust)

dev.off()

## 5. Affinity propagation (AP), too slow to be realized
cat(".......... No.5 Affinity propagation (AP), too slow to be realized ..........\n", sep="")

## 6. Gap statistics
cat(".......... No.6 Gap Statistics for Estimation ..........\n", sep="")
clusGap(TOM, kmeans, 15, B = 100, verbose = interactive())
## graphics still needed

## 7. Clustergrams to visualize cluster assignment

## 8. NbClust pacakges
cat(".......... No.8 NbClust ..........\n", sep="")
figureFile = paste(FigureDir, "03.Kmeans.TOM.ClusterNo.NbClust.emf", sep="")
emf(file = figureFile, 
	width = 6, 
	height = 5, 
	family="Helvetica"
)

nb = NbClust(TOM, 
		diss = "NULL", 
		distance = "euclidean", 
		min.nc = 2, 
		max.nc = 15, 
		method = "kmeans", 
		index = "alllong", 
		alphaBeale = 0.1)
hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))

dev.off()

## 9. for high-dimension data
cat(".......... No.9 High-dimension data with boostrap resampling ..........\n", sep="")

figureFile = paste(FigureDir, "03.Kmeans.TOM.ClusterNo.Pv.emf", sep="")
emf(file = figureFile, 
	width = 6, 
	height = 5, 
	family="Helvetica"
)

pv = pvclust(TOM)
plot(pv)

dev.off()

rm(list = ls())
gc()

