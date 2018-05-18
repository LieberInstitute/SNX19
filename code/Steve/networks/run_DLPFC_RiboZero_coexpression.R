#qsub -l bluejay -l mf=45G,h_vmem=55G,h_stack=256M -M stephensemick@gmail.com -cwd -b y R CMD BATCH --no-save run_DLPFC_RiboZero_coexpression.R
### 
library(GenomicRanges)
library(limma)
library(matrixStats)
library(RColorBrewer)
library(jaffelab)
library(jaffelab)
library(SummarizedExperiment)
###

## functions
getRPKM = recount::getRPKM
getRPM = function(rse, target = 10e6) {
	bg = matrix(rep(colSums(assays(rse)$counts)/target),
		nc = ncol(rse), nr = nrow(rse),byrow=TRUE)
	assays(rse)$counts/bg
}

## load counts
load("/dcl01/lieber/RNAseq/Datasets/DLPFC_Ribozero/rawCounts_DLPFC_Ribozero.rda")
theChrs = paste0("chr", c(1:22,"X","Y","MT"))

## make expression sets
pdDF = DataFrame(pd)

## gene
geneMapGR = makeGRangesFromDataFrame(geneMap, keep=TRUE)
rse_gene = SummarizedExperiment(
	assays = list('counts' = as.matrix(geneCounts)),
    colData = pdDF, rowRanges = geneMapGR)

## gene
exonMapGR = makeGRangesFromDataFrame(exonMap, keep=TRUE)
rse_exon = SummarizedExperiment(
	assays = list('counts' = as.matrix(exonCounts)),
    colData = pdDF, rowRanges = exonMapGR)
	
## junction
jCounts = as.matrix(as.data.frame(jCounts))
jIndex = which(rowSums(jCounts > 0) > 2)
jCounts = jCounts[jIndex,]
jMap = jMap[jIndex]
rse_jxn = SummarizedExperiment(
	assays = list('counts' = jCounts),
    colData = pdDF, rowRanges = jMap)
############
## filter samples to postnatal
keepIndex = which(colData(rse_gene)$Race %in% c("CAUC", "AA") & 
	colData(rse_gene)$age > 13)
rse_gene = rse_gene[,keepIndex]
rse_exon = rse_exon[,keepIndex]
rse_jxn = rse_jxn[,keepIndex]
# filter
geneIndex = rowMeans(getRPKM(rse_gene, length = "Length")) > 0.105
rse_gene = rse_gene[geneIndex,]
yGene = as.matrix(log2(geneRpkm+1))

exonIndex =  rowMeans(getRPKM(rse_exon, length = "Length")) > 0.1 &
				seqnames(rse_exon) %in% theChrs 
rse_exon = rse_exon[exonIndex,]
yExon = as.matrix(log2(exonRpkm+1))

jxnIndex =  rowMeans(getRPM(rse_jxn)) > 0.5 & rowData(rse_jxn)$code != "Novel" 
rse_jxn = rse_jxn[jxnIndex,]
## joint annotation
geneMap$Class = "InEns"
geneMap$meanExprs = rowMeans(geneRpkm)
geneMap$Geneid = rownames(geneMap)
exonMap$Class = "InEns"
exonMap$meanExprs = rowMeans(exonRpkm)

tmp = as.data.frame(jMap)
colnames(tmp)[1:3] = c("Chr", "Start", "End")
jMapDf = tmp[,c("Chr", "Start", "End", "code", 
	"newGeneID", "newGeneSymbol")]
jMapDf$meanExprs= rowMeans(jRpkm)
colnames(jMapDf)[4:6] = c("Class", "Geneid", "Symbol")

map = rbind(geneMap[,colnames(jMapDf)],
	exonMap[,colnames(jMapDf)], jMapDf)
map$featureID = rownames(map)
map$Type = rep(c("Gene", "Exon", "Junction"),
	c(nrow(geneMap), nrow(exonMap), nrow(jMapDf)))
###################
## bind and PCA
yExprs = rbind(yGene, yExon, yJxn)

# oo = order(rowSds(yExprs),decreasing=TRUE)[1:50000]
# pca = prcomp(t(yExprs[oo,]))

## add in qSVs
degCovAdj = read.csv("/users/ajaffe/Lieber/Projects/RNAseq/Consortium/PhaseI/degradationMat_DLPFC_polyA_Phase1.csv.gz",
	as.is=TRUE,row.names=1)
degCovAdj = degCovAdj[,colnames(yExprs)]
qSVs = prcomp(t(log2(degCovAdj+1)))$x[,1:15]

######################
# snx19 d9 expression
features = read.csv('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest.csv', stringsAsFactors=FALSE)
theJuncs = features$Feature
names(theJuncs) <- features$Tag
theJuncs = theJuncs[theJuncs%in% rownames(yExprs)] #dropping foi not available in this data
## filter for age
aIndex = which(pd$Age > 13)
yExprs2 = yExprs[,aIndex]
degCovAdj2 = degCovAdj[,aIndex]
pd2 = pd[aIndex,]
	
gIndexes = splitit(paste0(pd2$Dx, "_", pd2$Race)) #eQTL analysis for four subgroups: Caucasian controls, African American controls, Caucasian Schizophrenics, and African American Schizophrenics.
gIndexes$CAUC = which(pd2$Race == "CAUC") #eQTL analysis for Caucasians
gIndexes$AA = which(pd2$Race == "AA") #eQTL analysis fo   r African Americans
gIndexes$All =1:nrow(pd2)
tIndexes = splitit(map$Type)

library(limma)	
library(parallel)
statList = mclapply(gIndexes, function(ii) { 
	tmpList = lapply(theJuncs, function(jj) {
		cat(".")
		yy = yExprs2[jj,ii]
		qSVs2 = prcomp(t(log2(degCovAdj[,ii]+1)))$x[,1:10]
		mod2 = model.matrix(~as.numeric(yy) +
			pd2$Sex[ii] + pd2$Age[ii] +qSVs2)
		f = lmFit(yExprs2[,ii], mod2)
		e = ebayes(f)
		out = data.frame(log2FC = f$coef[,2], 
			tstat = e$t[,2], pval = e$p[,2])
		## add FDR
		out$fdr = NA
		for(k in seq(along=tIndexes)) {
			kk = tIndexes[[k]]
			out$fdr[kk] = p.adjust(out$pval[kk],"fdr")
		}
		return(out)
	})
	stats = do.call("cbind", tmpList)
},mc.cores=1)
outStats = do.call("cbind", statList)
outStats = cbind(map,outStats)

save(outStats, file = "rdas/Multivariate_Coexpression_SNX19_RiboZero_FullStats.rda" )
write.csv(outStats, file=gzfile("tables/Multivariate_Coexpression_SNX19_RiboZero_FullStats.csv.gz"),
	row.names=FALSE)
table(rowSums(outStats[,grep("pval", colnames(outStats))] < 0.01) )
