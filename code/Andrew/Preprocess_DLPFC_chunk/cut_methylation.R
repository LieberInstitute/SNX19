###
###
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

library(MatrixEQTL)
library(sva)
library(minfi)
library(GenomicRanges)

## load data
load("/dcs01/ajaffe/Brain/DNAm/ECD2014/rdas/meQTL_data.rda")
colnames(p) = rownames(pd)

## take 2 Mb region
gr = GRanges("chr11", IRanges(129718630, 131718630))
oo = findOverlaps(gr, map)

methMat = p[subjectHits(oo),]

methMap = as.data.frame(map[subjectHits(oo)])
colnames(methMap)[1:2] = c("chr","pos")
methMap = methMap[,-(3:5)]
methPheno = pd
save(methMat, methMap, methPheno, file="SNX19_DNAm_data_DLPFC.rda")