
##########################
#Author: Liang Ma, Qiang Chen
#SNX19
#calculating 9 groups
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
mod = model.matrix(~pd2$snpPC1 + pd2$snpPC2 + pd2$snpPC3 + pd2$snpPC4 + pd2$snpPC5 + exonPCs)
colnames(mod)[2:ncol(mod)] = c(paste0("snpPC",1:5), paste0("exprsPC",1:15))

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

##################
# all
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


# all CAUC 
#-------------------------------
caucIndex=which(pd2$Race=="CAUC")

snpSliceCauc = SlicedData$new(snp2[,caucIndex])
snpSliceCauc$ResliceCombined(sliceSize = 5000)

exprsCauc = SlicedData$new(rbind(geneRpkm2[,caucIndex], exonRpkm2[,caucIndex], jRp80M2[,caucIndex]))

modCauc = mod[caucIndex,]
covsCauc = SlicedData$new(t(modCauc[,-1]))

meJointCauc = Matrix_eQTL_main(snps=snpSliceCauc, gene = exprsCauc, 
	cvrt = covsCauc, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = pos, 
	useModel = modelLINEAR,	cisDist=5e6,
	pvalue.hist = 100,min.pv.by.genesnp = T)	

# all aa 
#----------------------------
aaIndex=which(pd2$Race=="AA")

snpSliceaa = SlicedData$new(snp2[,aaIndex])
snpSliceaa$ResliceCombined(sliceSize = 5000)

exprsaa = SlicedData$new(rbind(geneRpkm2[,aaIndex], exonRpkm2[,aaIndex], jRp80M2[,aaIndex]))

modaa = mod[aaIndex,]
covsaa = SlicedData$new(t(modaa[,-1]))

meJoint.aa = Matrix_eQTL_main(snps=snpSliceaa, gene = exprsaa, 
                              cvrt = covsaa, output_file_name.cis =  ".ctxt" ,
                              pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                              snpspos = snpspos, genepos = pos, 
                              useModel = modelLINEAR,  cisDist=5e6,
                              pvalue.hist = 100,min.pv.by.genesnp = T)  
# all ct 
#----------------------------
ctIndex=which(pd2$Dx=="Control")

modct = mod[ctIndex,]
covsct = SlicedData$new(t(modct[,-1]))

snpSlicect = SlicedData$new(snp2[,ctIndex])
snpSlicect$ResliceCombined(sliceSize = 5000)

exprsct = SlicedData$new(rbind(geneRpkm2[,ctIndex], exonRpkm2[,ctIndex], jRp80M2[,ctIndex]))

meJoint.ct = Matrix_eQTL_main(snps=snpSlicect, gene = exprsct, 
                              cvrt = covsct, output_file_name.cis =  ".ctxt" ,
                              pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                              snpspos = snpspos, genepos = pos, 
                              useModel = modelLINEAR,  cisDist=5e6,
                              pvalue.hist = 100,min.pv.by.genesnp = T)	

# all sz
#----------------------------
szIndex=which(pd2$Dx=="Schizo")

modsz = mod[szIndex,]
covssz = SlicedData$new(t(modsz[,-1]))

snpSlicesz = SlicedData$new(snp2[,szIndex])
snpSlicesz$ResliceCombined(sliceSize = 5000)

exprssz = SlicedData$new(rbind(geneRpkm2[,szIndex], exonRpkm2[,szIndex], jRp80M2[,szIndex]))

meJoint.sz = Matrix_eQTL_main(snps=snpSlicesz, gene = exprssz, 
                              cvrt = covssz, output_file_name.cis =  ".ctxt" ,
                              pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                              snpspos = snpspos, genepos = pos, 
                              useModel = modelLINEAR,  cisDist=5e6,
                              pvalue.hist = 100,min.pv.by.genesnp = T)	
# all cauc.ct
#----------------------------
cauc.ctIndex=which((pd2$Race=="CAUC") & (pd2$Dx=="Control"))

modcauc.ct = mod[cauc.ctIndex,]
covscauc.ct = SlicedData$new(t(modcauc.ct[,-1]))

snpSlicecauc.ct = SlicedData$new(snp2[,cauc.ctIndex])
snpSlicecauc.ct$ResliceCombined(sliceSize = 5000)

exprscauc.ct = SlicedData$new(rbind(geneRpkm2[,cauc.ctIndex], exonRpkm2[,cauc.ctIndex], jRp80M2[,cauc.ctIndex]))

meJoint.cauc.ct = Matrix_eQTL_main(snps=snpSlicecauc.ct, gene = exprscauc.ct, 
                                   cvrt = covscauc.ct, output_file_name.cis =  ".ctxt" ,
                                   pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                                   snpspos = snpspos, genepos = pos, 
                                   useModel = modelLINEAR,  cisDist=5e6,
                                   pvalue.hist = 100,min.pv.by.genesnp = T)	
# all cauc.sz
#----------------------------
cauc.szIndex=which((pd2$Race=="CAUC") & (pd2$Dx=="Schizo"))

modcauc.sz = mod[cauc.szIndex,]
covscauc.sz = SlicedData$new(t(modcauc.sz[,-1]))

snpSlicecauc.sz = SlicedData$new(snp2[,cauc.szIndex])
snpSlicecauc.sz$ResliceCombined(sliceSize = 5000)

exprscauc.sz = SlicedData$new(rbind(geneRpkm2[,cauc.szIndex], exonRpkm2[,cauc.szIndex], jRp80M2[,cauc.szIndex]))

meJoint.cauc.sz = Matrix_eQTL_main(snps=snpSlicecauc.sz, gene = exprscauc.sz, 
                                   cvrt = covscauc.sz, output_file_name.cis =  ".ctxt" ,
                                   pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                                   snpspos = snpspos, genepos = pos, 
                                   useModel = modelLINEAR,  cisDist=5e6,
                                   pvalue.hist = 100,min.pv.by.genesnp = T)	
# all aa.ct
#----------------------------
aa.ctIndex=which((pd2$Race=="AA") & (pd2$Dx=="Control"))

modaa.ct = mod[aa.ctIndex,]
covsaa.ct = SlicedData$new(t(modaa.ct[,-1]))

snpSliceaa.ct = SlicedData$new(snp2[,aa.ctIndex])
snpSliceaa.ct$ResliceCombined(sliceSize = 5000)

exprsaa.ct = SlicedData$new(rbind(geneRpkm2[,aa.ctIndex], exonRpkm2[,aa.ctIndex], jRp80M2[,aa.ctIndex]))

meJoint.aa.ct = Matrix_eQTL_main(snps=snpSliceaa.ct, gene = exprsaa.ct, 
                                 cvrt = covsaa.ct, output_file_name.cis =  ".ctxt" ,
                                 pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                                 snpspos = snpspos, genepos = pos, 
                                 useModel = modelLINEAR,  cisDist=5e6,
                                 pvalue.hist = 100,min.pv.by.genesnp = T)  
# all aa.sz
#----------------------------
aa.szIndex=which((pd2$Race=="AA") & (pd2$Dx=="Schizo"))


modaa.sz = mod[aa.szIndex,]
covsaa.sz = SlicedData$new(t(modaa.sz[,-1]))

snpSliceaa.sz = SlicedData$new(snp2[,aa.szIndex])
snpSliceaa.sz$ResliceCombined(sliceSize = 5000)

exprsaa.sz = SlicedData$new(rbind(geneRpkm2[,aa.szIndex], exonRpkm2[,aa.szIndex], jRp80M2[,aa.szIndex]))

meJoint.aa.sz = Matrix_eQTL_main(snps=snpSliceaa.sz, gene = exprsaa.sz, 
                                 cvrt = covsaa.sz, output_file_name.cis =  ".ctxt" ,
                                 pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
                                 snpspos = snpspos, genepos = pos, 
                                 useModel = modelLINEAR,  cisDist=5e6,
                                 pvalue.hist = 100,min.pv.by.genesnp = T)  


#save
#-------------------
save()




