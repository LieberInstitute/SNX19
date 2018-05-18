## Create correlation matrix for each dataset for 10 features of interest
#load features of interest
load('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild_03_22_18.rda')

#case_control_differential_expression_5_regions.R
library(SummarizedExperiment)
library(limma)

####====== DLPFC PolyA ======#### 
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_polyA/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_n495_updateMap_03082017.rda")

## Filtering people
keepIndex = which(pd$Race %in% c("CAUC", "AA") &  pd$age > 13)
pd=pd[keepIndex, ]
geneRpkm = geneRpkm[,keepIndex]
exonRpkm = exonRpkm[,keepIndex]
jRpkm = jRpkm[,keepIndex]

#expression object
gIndex=which(rowMeans(geneRpkm) > 0.01) #selecting out genes with sufficient expression (genes with less than 0.01 read per kilobase million are not perserved)
geneRpkm = geneRpkm[gIndex,]
geneMap = geneMap[gIndex,]

eIndex=which(rowMeans(exonRpkm) > 0.01)  #selecting out exons with sufficient expression (exons with less than 0.01 read per kilobase million are not perserved)
exonRpkm = exonRpkm[eIndex,]
exonMap = exonMap[eIndex, ]

jIndex = which(rowMeans(jRpkm) > 0.2)#There are an excessive number of junctions. Apply a stricter filter and drop novel junctions to reduce the number.
jRpkm = jRpkm[jIndex,]
jMap = jMap[jIndex]

yExprs_DlpfcPolya = log2(rbind(geneRpkm, exonRpkm, jRpkm)+1)
pd_DlpfcPolyA = pd

####====== DLPFC RiboZero ======#### 
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_RiboZero/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_RiboZero_n485_updateMap.rda")

## Filtering people
keepIndex = which(pd$Race %in% c("CAUC", "AA") &  pd$age > 13)
pd=pd[keepIndex, ]
geneRpkm = geneRpkm[,keepIndex]
exonRpkm = exonRpkm[,keepIndex]
jRpkm = jRpkm[,keepIndex]

#expression object
gIndex=which(rowMeans(geneRpkm) > 0.01) #selecting out genes with sufficient expression (genes with less than 0.01 read per kilobase million are not perserved)
geneRpkm = geneRpkm[gIndex,]
geneMap = geneMap[gIndex,]

eIndex=which(rowMeans(exonRpkm) > 0.01)  #selecting out exons with sufficient expression (exons with less than 0.01 read per kilobase million are not perserved)
exonRpkm = exonRpkm[eIndex,]
exonMap = exonMap[eIndex, ]

jIndex = which(rowMeans(jRpkm) > 0.2)#There are an excessive number of junctions. Apply a stricter filter and drop novel junctions to reduce the number.
jRpkm = jRpkm[jIndex,]
jMap = jMap[jIndex]


yExprs_DlpfcRiboZero = log2(rbind(geneRpkm, exonRpkm, jRpkm)+1)
pd_DlpfcRiboZero = pd

####====== Caudate ======#### 
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/Caudate_RiboZero/rawAndRpkmCounts_plusGenotype_SNX19_Caudate_n468_updateMap.rda")

## Filtering people
keepIndex = which(pd$Race %in% c("CAUC", "AA") &  pd$age > 13)
pd=pd[keepIndex, ]
geneRpkm = geneRpkm[,keepIndex]
exonRpkm = exonRpkm[,keepIndex]
jRpkm = jRpkm[,keepIndex]

#expression object
gIndex=which(rowMeans(geneRpkm) > 0.01) #selecting out genes with sufficient expression (genes with less than 0.01 read per kilobase million are not perserved)
geneRpkm = geneRpkm[gIndex,]
geneMap = geneMap[gIndex,]

eIndex=which(rowMeans(exonRpkm) > 0.01)  #selecting out exons with sufficient expression (exons with less than 0.01 read per kilobase million are not perserved)
exonRpkm = exonRpkm[eIndex,]
exonMap = exonMap[eIndex, ]

jIndex = which(rowMeans(jRpkm) > 0.2)#There are an excessive number of junctions. Apply a stricter filter and drop novel junctions to reduce the number.
jRpkm = jRpkm[jIndex,]
jMap = jMap[jIndex]

yExprs_Caudate = log2(rbind(geneRpkm, exonRpkm, jRpkm)+1)
pd_Caudate = pd

####====== Hippo ======#### 
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/HIPPO_RiboZero/rawAndRpkmCounts_plusGenotype_SNX19_HIPPO_RiboZero_n452_updateMap.rda")

## Filtering people
keepIndex = which(pd$Race %in% c("CAUC", "AA") &  pd$Age > 13)
pd=pd[keepIndex, ]
geneRpkm = geneRpkm[,keepIndex]
exonRpkm = exonRpkm[,keepIndex]
jRpkm = jRpkm[,keepIndex]

#expression object
gIndex=which(rowMeans(geneRpkm) > 0.01) #selecting out genes with sufficient expression (genes with less than 0.01 read per kilobase million are not perserved)
geneRpkm = geneRpkm[gIndex,]
geneMap = geneMap[gIndex,]

eIndex=which(rowMeans(exonRpkm) > 0.01)  #selecting out exons with sufficient expression (exons with less than 0.01 read per kilobase million are not perserved)
exonRpkm = exonRpkm[eIndex,]
exonMap = exonMap[eIndex, ]

jIndex = which(rowMeans(jRpkm) > 0.2)#There are an excessive number of junctions. Apply a stricter filter and drop novel junctions to reduce the number.
jRpkm = jRpkm[jIndex,]
jMap = jMap[jIndex]


yExprs_Hippo = log2(rbind(geneRpkm, exonRpkm, jRpkm)+1)
pd_Hippo = pd

####====== CMC DLPFC ======#### 
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/CMC/rawAndRpkmCounts_plusGenotype_SNX19_CMC_n547.rda")

## filter samples to postnatal
pd$Age_of_Death[pd$Age_of_Death=="90+"] <- NA
pd <- plyr::rename(pd, c("Age_of_Death"="age"))
pd$age <- as.numeric(pd$age)
pd$Race = plyr::mapvalues(pd$Ethnicity, c("African-American", "Caucasian"), c("AA","CAUC") )
pd$Dx = plyr::mapvalues(pd$Dx, c("SCZ"), c("Schizo") )
pd$Sex = plyr::mapvalues(pd$Gender, c("Male", "Female"), c("M","F") )

## Filtering people
keepIndex = which(pd$Race %in% c("CAUC", "AA") &  pd$age > 13)
pd=pd[keepIndex, ]
geneRpkm = geneRpkm[,keepIndex]
exonRpkm = exonRpkm[,keepIndex]
jRpkm = jRpkm[,keepIndex]

#expression object
gIndex=which(rowMeans(geneRpkm) > 0.01) #selecting out genes with sufficient expression (genes with less than 0.01 read per kilobase million are not perserved)
geneRpkm = geneRpkm[gIndex,]
geneMap = geneMap[gIndex,]

eIndex=which(rowMeans(exonRpkm) > 0.01)  #selecting out exons with sufficient expression (exons with less than 0.01 read per kilobase million are not perserved)
exonRpkm = exonRpkm[eIndex,]
exonMap = exonMap[eIndex, ]

jIndex = which(rowMeans(jRpkm) > 0.2)#There are an excessive number of junctions. Apply a stricter filter and drop novel junctions to reduce the number.
jRpkm = jRpkm[jIndex,]
jMap = jMap[jIndex]

yExprs_CMC = log2(rbind(geneRpkm, exonRpkm, jRpkm)+1)
pd_CMC = pd

####### calculate correlation for each dataset
library(tibble)
library(Hmisc)

## Dlpfc PolyA
corFOI_DlpfcPolyA <- rcorr(t(yExprs_DlpfcPolya[foi$Feature,]), type="pearson")

corFOI_DlpfcPolyA_r = corFOI_DlpfcPolyA$r
corFOI_DlpfcPolyA_r[upper.tri(corFOI_DlpfcPolyA_r, diag=F )] <-NA
rownames(corFOI_DlpfcPolyA_r) = colnames(corFOI_DlpfcPolyA_r) = foi[match(rownames(corFOI_DlpfcPolyA_r),foi$Feature ),'Tag']
corFOI_DlpfcPolyA_r = rownames_to_column(as.data.frame(corFOI_DlpfcPolyA_r),var='FeatureName')

corFOI_DlpfcPolyA_p = corFOI_DlpfcPolyA$P
corFOI_DlpfcPolyA_p[upper.tri(corFOI_DlpfcPolyA_p, diag=F )] <-NA
rownames(corFOI_DlpfcPolyA_p) = colnames(corFOI_DlpfcPolyA_p) = foi[match(rownames(corFOI_DlpfcPolyA_p),foi$Feature ),'Tag']
corFOI_DlpfcPolyA_p = rownames_to_column(as.data.frame(corFOI_DlpfcPolyA_p),var='FeatureName')

## Dlpfc RiboZero
corFOI_DlpfcRiboZero= rcorr(t(yExprs_DlpfcRiboZero[foi$Feature,]), type="pearson")

corFOI_DlpfcRiboZero_r = corFOI_DlpfcRiboZero$r
corFOI_DlpfcRiboZero_r[upper.tri(corFOI_DlpfcRiboZero_r, diag=F )] <-NA
rownames(corFOI_DlpfcRiboZero_r) = colnames(corFOI_DlpfcRiboZero_r) = foi[match(rownames(corFOI_DlpfcRiboZero_r),foi$Feature ),'Tag']
corFOI_DlpfcRiboZero_r = rownames_to_column(as.data.frame(corFOI_DlpfcRiboZero_r),var='FeatureName')

corFOI_DlpfcRiboZero_p = corFOI_DlpfcRiboZero$P
corFOI_DlpfcRiboZero_p[upper.tri(corFOI_DlpfcRiboZero_p, diag=F )] <-NA
rownames(corFOI_DlpfcRiboZero_p) = colnames(corFOI_DlpfcRiboZero_p) = foi[match(rownames(corFOI_DlpfcRiboZero_p),foi$Feature ),'Tag']
corFOI_DlpfcRiboZero_p = rownames_to_column(as.data.frame(corFOI_DlpfcRiboZero_p),var='FeatureName')

## CMC
corFOI_CMC= rcorr(t(yExprs_CMC[foi$Feature,]), type="pearson")

corFOI_CMC_r = corFOI_CMC$r
corFOI_CMC_r[upper.tri(corFOI_CMC_r, diag=F )] <-NA
rownames(corFOI_CMC_r) = colnames(corFOI_CMC_r) = foi[match(rownames(corFOI_CMC_r),foi$Feature ),'Tag']
corFOI_CMC_r = rownames_to_column(as.data.frame(corFOI_CMC_r),var='FeatureName')

corFOI_CMC_p = corFOI_CMC$P
corFOI_CMC_p[upper.tri(corFOI_CMC_p, diag=F )] <-NA
rownames(corFOI_CMC_p) = colnames(corFOI_CMC_p) = foi[match(rownames(corFOI_CMC_p),foi$Feature ),'Tag']
corFOI_CMC_p = rownames_to_column(as.data.frame(corFOI_CMC_p),var='FeatureName')

## Hippo
corFOI_Hippo= rcorr(t(yExprs_Hippo[foi$Feature,]), type ='pearson')

corFOI_Hippo_r = corFOI_Hippo$r
corFOI_Hippo_r[upper.tri(corFOI_Hippo_r, diag=F )] <-NA
rownames(corFOI_Hippo_r) = colnames(corFOI_Hippo_r) = foi[match(rownames(corFOI_Hippo_r),foi$Feature ),'Tag']
corFOI_Hippo_r = rownames_to_column(as.data.frame(corFOI_Hippo_r),var='FeatureName')

corFOI_Hippo_p = corFOI_Hippo$P
corFOI_Hippo_p[upper.tri(corFOI_Hippo_p, diag=F )] <-NA
rownames(corFOI_Hippo_p) = colnames(corFOI_Hippo_p) = foi[match(rownames(corFOI_Hippo_p),foi$Feature ),'Tag']
corFOI_Hippo_p = rownames_to_column(as.data.frame(corFOI_Hippo_p),var='FeatureName')


## Caudate
corFOI_Caudate= rcorr(t(yExprs_Caudate[foi$Feature,]), type ='pearson')

corFOI_Caudate_r = corFOI_Caudate$r
corFOI_Caudate_r[upper.tri(corFOI_Caudate_r, diag=F )] <-NA
rownames(corFOI_Caudate_r) = colnames(corFOI_Caudate_r) = foi[match(rownames(corFOI_Caudate_r),foi$Feature ),'Tag']
corFOI_Caudate_r = rownames_to_column(as.data.frame(corFOI_Caudate_r),var='FeatureName')

corFOI_Caudate_p = corFOI_Caudate$P
corFOI_Caudate_p[upper.tri(corFOI_Caudate_p, diag=F )] <-NA
rownames(corFOI_Caudate_p) = colnames(corFOI_Caudate_p) = foi[match(rownames(corFOI_Caudate_p),foi$Feature ),'Tag']
corFOI_Caudate_p = rownames_to_column(as.data.frame(corFOI_Caudate_p),var='FeatureName')


openxlsx::write.xlsx(list(Dlpfc_LIBD_PolyA=corFOI_DlpfcPolyA_r,
						   Dlpfc_LIBD_RiboZero=corFOI_DlpfcRiboZero_r,
						   Dlpfc_CMC_RiboZero=corFOI_CMC_r,
						   Hippo_LIBD_RiboZero=corFOI_Hippo_r,
						   Caudate_LIBD_RiboZero=corFOI_Caudate_r), file='/dcl01/lieber/ajaffe/Steve/SNX19/tables/correlation_matrices_withinDataset_betweenFeature_pearson_R.xlsx')
openxlsx::write.xlsx(list(Dlpfc_LIBD_PolyA=corFOI_DlpfcPolyA_p,
						   Dlpfc_LIBD_RiboZero=corFOI_DlpfcRiboZero_p,
						   Dlpfc_CMC_RiboZero=corFOI_CMC_p,
						   Hippo_LIBD_RiboZero=corFOI_Hippo_p,
						   Caudate_LIBD_RiboZero=corFOI_Caudate_p), file='/dcl01/lieber/ajaffe/Steve/SNX19/tables/correlation_matrices_withinDataset_betweenFeature_pearson_P.xlsx')						   
						   
