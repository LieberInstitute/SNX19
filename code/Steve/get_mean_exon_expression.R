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
exonMap_DlpfcPolya=exonMap

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
exonMap_DlpfcRiboZero=exonMap

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
exonMap_Caudate=exonMap

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
exonMap_Hippo=exonMap
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
exonMap_CMC=exonMap

####### calculate correlation for each dataset
library(tibble)
library(Hmisc)

exonMap_All = rbind(exonMap_DlpfcPolya,exonMap_DlpfcRiboZero,exonMap_Caudate,exonMap_Hippo,exonMap_CMC)
exonMap_All$Coordinates = paste0(exonMap_All$Chr,":",exonMap_All$Start,"-",exonMap_All$End)
exonMap_All=exonMap_All[!duplicated(exonMap_All$Coordinates),]
exonMap_SNX19 = exonMap_All[exonMap_All$Symbol=="SNX19", ]

meanExprs = sapply(list(yExprs_DlpfcPolya,yExprs_DlpfcRiboZero,yExprs_Caudate,yExprs_Hippo,yExprs_CMC),function(x) rowMeans(x[rownames(exonMap_SNX19),]) )
colnames(meanExprs) <-  paste0("logScale_meanExprs_", c("DlpfcPolya","yExprs_DlpfcRiboZero","yExprs_Caudate","yExprs_Hippo","yExprs_CMC") )
meanExprs=as.data.frame(meanExprs)
meanExprs=cbind(exonMap_SNX19,meanExprs)

#### add ensembl exon ids
ensemblGTF = rtracklayer::import(con="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ensembl/Homo_sapiens.GRCh37.75.gtf", format="gtf")
ensEXONS = as.data.frame(ensemblGTF)[which(ensemblGTF$type=="exon"),c("seqnames","start","end","strand","gene_id","exon_id","transcript_id")]
names(ensEXONS) = c("Chr","Start","End","Strand","ensemblID","exon_ensemblID","transcript_id")
rownames(ensEXONS) = NULL
ensEXONS$eID = paste0("e", rownames(ensEXONS))
ensEXONS$SNX19_ID = paste0("chr", ensEXONS$Chr,":",ensEXONS$Start,"-", ensEXONS$End, "(",ensEXONS$Strand, ")")
meanExprs$txMatchID = rownames(meanExprs)
eeMatch = match(meanExprs$txMatchID,ensEXONS$SNX19_ID)
meanExprs$ensembl_exon_id <- ensEXONS[eeMatch,'exon_ensemblID']
meanExprs$ensembl_transcript_id <- ensEXONS[eeMatch,'transcript_id']

## add tag names
tagNames = read.csv('/dcl01/lieber/ajaffe/Steve/SNX19/snx19_exon_junction_tagnames.csv',stringsAsFactors=FALSE)
tagNames$TagName = tagNames$Feature
meanExprs$tagName = tagNames[match(rownames(meanExprs), tagNames$Position_hg19_strand ),'TagName']

absAdd= 2^meanExprs[, grep("logScale_meanExprs_",colnames(meanExprs) ) ]+1
colnames(absAdd) <- gsub("logScale","absScale",colnames(absAdd) )
meanExprs = cbind(meanExprs, absAdd )


write.csv(meanExprs,file='/dcl01/lieber/ajaffe/Steve/SNX19/tables/SNX19_mean_exon_expression_logRPKM.csv',row.names=FALSE)


###########
load('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild_03_22_18.rda')
juncPull = tagNames[grepl("junc",tagNames$Feature),c('Position_hg19_strand','TagName')]
colnames(juncPull) <- c("Feature","Tag")
addFoi = foi[foi$exprsType=="Junction",c('Feature','Tag')]
juncPull=rbind(juncPull,addFoi)
juncPull = juncPull[!duplicated(juncPull),]
juncPull = cbind(juncPull, as.data.frame(jMap)[match(juncPull$Feature, names(jMap)),] )


meanExprsJxn = sapply(list(yExprs_DlpfcPolya,yExprs_DlpfcRiboZero,yExprs_Caudate,yExprs_Hippo,yExprs_CMC),function(x) rowMeans(x[juncPull$Feature,]) )
colnames(meanExprsJxn) <- paste0("logScale_meanExprs_",c("DlpfcPolya","yExprs_DlpfcRiboZero","yExprs_Caudate","yExprs_Hippo","yExprs_CMC") )
meanExprsJxn=as.data.frame(meanExprsJxn)
meanExprsJxn=cbind(juncPull,meanExprsJxn)
absAdd= 2^meanExprsJxn[, grep("logScale_meanExprs_",colnames(meanExprsJxn) ) ]+1
colnames(absAdd) <- gsub("logScale","absScale",colnames(absAdd) )
meanExprsJxn = cbind(meanExprsJxn, absAdd )
meanExprsJxn$ensemblTx = sapply(meanExprsJxn$ensemblTx , function(x) paste(unlist(x),collapse=';') )
write.csv(meanExprsJxn,file='/dcl01/lieber/ajaffe/Steve/SNX19/tables/SNX19_mean_jxn_expression_logRPKM.csv',row.names=FALSE)
