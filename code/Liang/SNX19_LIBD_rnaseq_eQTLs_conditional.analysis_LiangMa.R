
##########################
#Author: Liang Ma, Qiang Chen
#SNX19
#conditional analysis for junc8.10 and junc8.8a
##########################

library(MatrixEQTL)
library(GenomicRanges)
library(sva)
library(XVector)
library(Biostrings)

#### load data
load("rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_n495.Wholegenome_TotalCounts.rda")

# make RPKM
#-------------------------
bgGene = matrix(rep(pd$totalMapped), nc = nrow(pd), nr = nrow(geneCounts),	byrow=TRUE)
widGGene = matrix(rep(geneMap$Length), nr = nrow(geneCounts), nc = nrow(pd),	byrow=FALSE)
geneRpkm = geneCounts/(widGGene/1000)/(bgGene/1e6)


bgJunc = matrix(rep(pd$totalMapped/80e6), nc = nrow(pd), nr = nrow(jCounts), byrow=TRUE)
jRp80M = jCounts/bgJunc

## filter expression data
#-------------------------------
gIndex=which(rowMeans(geneRpkm) > 0.01)
geneRpkm = geneRpkm[gIndex,]
geneMap = geneMap[gIndex,]

eIndex=which(rowMeans(exonRpkm) > 0.01)
exonRpkm = exonRpkm[eIndex,]
exonMap = exonMap[eIndex,]

# jIndex = which(rowMeans(jRp80M) > 0.2 & jMap$code != "Novel")
# jRp80M = jRp80M[jIndex,]
# jMap= jMap[jIndex]

# filter 
#--------------------------
aIndex = which(pd$Age>13)

pd2 = pd[aIndex,]
snp2 = as.matrix(snpAll[,aIndex])
geneRpkm2 = as.matrix(log2(geneRpkm[,aIndex]+1))
exonRpkm2 = as.matrix(log2(exonRpkm[,aIndex]+1))
jRp80M2 = as.matrix(log2(jRp80M[,aIndex]+1))

exprs2 = rbind(geneRpkm2, exonRpkm2, jRp80M2);

# PCA
#-------------------
load("exonPCs_13plus.rda")



## model ####
#--------------------------

junc8.10 = exprs2[which(rownames(exprs2)=="chr11:130749607-130773149(*)"),]
junc8.8a = exprs2[which(rownames(exprs2)=="chr11:130765052-130773149(*)"),]

junc8.8a = data.frame(t(junc8.8a))
junc8.8a = as.matrix(junc8.8a)
junc8.10 = data.frame(t(junc8.10))
junc8.10 = as.matrix(junc8.10)


mod = model.matrix(~pd2$snpPC1 + pd2$snpPC2 + pd2$snpPC3 + pd2$snpPC4 + pd2$snpPC5 + exonPCs + junc8.8a)
colnames(mod)[2:ncol(mod)] = c(paste0("snpPC",1:5), paste0("exprsPC",1:15), paste0("junc8.8a"))

# # covariates
covs = SlicedData$new(t(mod[,-1])) 

###### snp position
snpspos = snpMapAll[,c("SNP","CHR","POS")] 
snpspos$CHR = paste0("chr",snpspos$CHR) 
colnames(snpspos) = c("name","chr","pos")

###### gene  position
posGene = geneMap[,1:3]
posGene$name = rownames(geneMap)
posGene = posGene[,c(4,1:3)]

##### exon  position 
posExon = exonMap[,2:4]
posExon$name = rownames(exonMap)
posExon = posExon[,c(4,1:3)]

##### junction  position 
posJxn = as.data.frame(jMap)[,1:3]
posJxn$name = names(jMap)
posJxn = posJxn[,c(4,1:3)]
names(posJxn)[2:4] = c("Chr", "Start","End")


####################
snpSlice = SlicedData$new(snp2)

snpSlice$ResliceCombined(sliceSize = 5000)

exprsSlice = SlicedData$new(rbind(geneRpkm2, exonRpkm2, jRp80M2))

pos = rbind(posGene, posExon, posJxn)

meJoint = Matrix_eQTL_main(snps=snpSlice, gene = exprsSlice, 
	cvrt = covs, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = pos, 
	useModel = modelLINEAR,	cisDist=5e6,
	pvalue.hist = 100, min.pv.by.genesnp = T)	


#save
#-------------------
save()






