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

load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/DLPFC_polyA/rdas/all3Feature_PCs_DLPFC_polyA_age13_matchedUp.rda')

yExprs_DlpfcPolya = rbind(geneRpkm, exonRpkm, jRpkm)

####### merge map objects

### common expression map to annotate
geneMap$ensemblID = rownames(geneMap)
geneMap$Type = 'Gene'
geneMap$Class = 'InEns'

exonMap$Class = "InEns"
exonMap$Type = "Exon"
exonMap = plyr::rename(exonMap, c('Geneid'='ensemblID') )


jMap=as.data.frame(jMap)
jMap$Type = "Junction"
jMap = plyr::rename(jMap, c('seqnames'='Chr','start'='Start','end'='End', 'strand'='Strand', 'newGeneID'='ensemblID', 'newGeneSymbol'='Symbol','code'='Class') )

name = c("ensemblID", "Symbol", "Type", "Class", 'Chr','Start','End', 'Strand')
exprsMap = rbind(geneMap[,name],
				 exonMap[,name],
				 jMap[,name])

#####				 
mod = model.matrix(~pd$Dx + pd$RIN + pd$Sex + pd$age + pd$snpPC1 + pd$snpPC2 + pd$snpPC3 + pd$snpPC4 + pd$snpPC5 +all3PCs[pd$RNum,1:15] ) #forming the model matrix from the SNP PCs and the expression PCs

fit = lmFit(yExprs_DlpfcPolya, mod)
fit <- eBayes(fit)
DlpfcPolyA_Stats = topTable(fit, coef=2, n=Inf, sort.by = 'none', genelist=exprsMap)
colnames(DlpfcPolyA_Stats) <- paste0("DlpfcPolyA_", colnames(DlpfcPolyA_Stats))

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

####### merge map objects

### common expression map to annotate
geneMap$ensemblID = rownames(geneMap)
geneMap$Type = 'Gene'
geneMap$Class = 'InEns'

exonMap$Class = "InEns"
exonMap$Type = "Exon"
exonMap = plyr::rename(exonMap, c('Geneid'='ensemblID') )


jMap=as.data.frame(jMap)
jMap$Type = "Junction"
jMap = plyr::rename(jMap, c('seqnames'='Chr','start'='Start','end'='End', 'strand'='Strand', 'newGeneID'='ensemblID', 'newGeneSymbol'='Symbol','code'='Class') )

name = c("ensemblID", "Symbol", "Type", "Class", 'Chr','Start','End', 'Strand')
exprsMap = rbind(geneMap[,name],
				 exonMap[,name],
				 jMap[,name])
	 
#####			
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/DLPFC_RiboZero/rdas/all3Feature_PCs_DLPFC_RiboZero_age13_matchedUp.rda')	 
mod = model.matrix(~pd$Dx + pd$RIN + pd$Sex + pd$age + pd$snpPC1 + pd$snpPC2 + pd$snpPC3 + pd$snpPC4 + pd$snpPC5 +all3PCs[pd$RNum,1:15] ) #forming the model matrix from the SNP PCs and the expression PCs

fit = lmFit(yExprs_DlpfcRiboZero, mod)
fit <- eBayes(fit)
DlpfcRiboZero_Stats = topTable(fit, coef=2, n=Inf, sort.by = 'none', genelist=exprsMap)
colnames(DlpfcRiboZero_Stats) <- paste0("DlpfcRiboZero_", colnames(DlpfcRiboZero_Stats))

table(DlpfcRiboZero_Stats$DlpfcRiboZero_adj.P.Val<0.05)

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

####### merge map objects

### common expression map to annotate
geneMap$ensemblID = rownames(geneMap)
geneMap$Type = 'Gene'
geneMap$Class = 'InEns'

exonMap$Class = "InEns"
exonMap$Type = "Exon"
exonMap = plyr::rename(exonMap, c('Geneid'='ensemblID') )


jMap=as.data.frame(jMap)
jMap$Type = "Junction"
jMap = plyr::rename(jMap, c('seqnames'='Chr','start'='Start','end'='End', 'strand'='Strand', 'newGeneID'='ensemblID', 'newGeneSymbol'='Symbol','code'='Class') )

name = c("ensemblID", "Symbol", "Type", "Class", 'Chr','Start','End', 'Strand')
exprsMap = rbind(geneMap[,name],
				 exonMap[,name],
				 jMap[,name])
	 
#####			
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Caudate/rdas/all3Feature_PCs_Caudate_RiboZero_age13_matchedUp.rda')	 
mod = model.matrix(~pd$Dx + pd$RIN + pd$Sex + pd$age + pd$snpPC1 + pd$snpPC2 + pd$snpPC3 + pd$snpPC4 + pd$snpPC5 +all3PCs[pd$RNum,1:15] ) #forming the model matrix from the SNP PCs and the expression PCs

fit = lmFit(yExprs_Caudate, mod)
fit <- eBayes(fit)
Caudate_Stats = topTable(fit, coef=2, n=Inf, sort.by = 'none', genelist=exprsMap)
colnames(Caudate_Stats) <- paste0("Caudate_", colnames(Caudate_Stats))

table(Caudate_Stats$Caudate_adj.P.Val<0.05)

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

####### merge map objects

### common expression map to annotate
geneMap$ensemblID = rownames(geneMap)
geneMap$Type = 'Gene'
geneMap$Class = 'InEns'

exonMap$Class = "InEns"
exonMap$Type = "Exon"
exonMap = plyr::rename(exonMap, c('Geneid'='ensemblID') )


jMap=as.data.frame(jMap)
jMap$Type = "Junction"
jMap = plyr::rename(jMap, c('seqnames'='Chr','start'='Start','end'='End', 'strand'='Strand', 'newGeneID'='ensemblID', 'newGeneSymbol'='Symbol','code'='Class') )

name = c("ensemblID", "Symbol", "Type", "Class", 'Chr','Start','End', 'Strand')
exprsMap = rbind(geneMap[,name],
				 exonMap[,name],
				 jMap[,name])
	 
#####			
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Hippo/rdas/all3Feature_PCs_Hippo_RiboZero_age13_matchedUp.rda')	 
mod = model.matrix(~pd$Dx + pd$RIN + pd$Sex + pd$Age + pd$snpPC1 + pd$snpPC2 + pd$snpPC3 + pd$snpPC4 + pd$snpPC5 +all3PCs[pd$RNum,1:15] ) #forming the model matrix from the SNP PCs and the expression PCs

fit = lmFit(yExprs_Hippo, mod)
fit <- eBayes(fit)
Hippo_Stats = topTable(fit, coef=2, n=Inf, sort.by = 'none', genelist=exprsMap)
colnames(Hippo_Stats) <- paste0("Hippo_", colnames(Hippo_Stats))

table(Hippo_Stats$Hippo_adj.P.Val<0.05)


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

####### merge map objects

### common expression map to annotate
geneMap$ensemblID = rownames(geneMap)
geneMap$Type = 'Gene'
geneMap$Class = 'InEns'

exonMap$Class = "InEns"
exonMap$Type = "Exon"
exonMap = plyr::rename(exonMap, c('Geneid'='ensemblID') )


jMap=as.data.frame(jMap)
jMap$Type = "Junction"
jMap = plyr::rename(jMap, c('seqnames'='Chr','start'='Start','end'='End', 'strand'='Strand', 'newGeneID'='ensemblID', 'newGeneSymbol'='Symbol','code'='Class') )

name = c("ensemblID", "Symbol", "Type", "Class", 'Chr','Start','End', 'Strand')
exprsMap = rbind(geneMap[,name],
				 exonMap[,name],
				 jMap[,name])
	 
#####			
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/CMC/rdas/all3Feature_PCs_CMC_DLPFC_age13_matchedUp.rda')	 
mod = model.matrix(~pd$Dx + pd$DLPFC_RNA_isolation_RIN + pd$Sex + pd$age + pd$snpPC1 + pd$snpPC2 + pd$snpPC3 + pd$snpPC4 + pd$snpPC5 +all3PCs[pd$DLPFC_RNA_Sequencing_Sample_ID,1:15] ) #forming the model matrix from the SNP PCs and the expression PCs

fit = lmFit(yExprs_CMC, mod)
fit <- eBayes(fit)
CMC_Stats = topTable(fit, coef=2, n=Inf, sort.by = 'none', genelist=exprsMap)
colnames(CMC_Stats) <- paste0("CMC_", colnames(CMC_Stats))

table(CMC_Stats$Hippo_adj.P.Val<0.05)

### Merge all results
ids = unique(c(rownames(DlpfcPolyA_Stats),
			 rownames(DlpfcRiboZero_Stats),
			 rownames(CMC_Stats),
			 rownames(Caudate_Stats),
			 rownames(Hippo_Stats) ) )
			 
mergedStats = cbind( DlpfcPolyA_Stats[match(ids,rownames(DlpfcPolyA_Stats)), ],
					 DlpfcRiboZero_Stats[match(ids,rownames(DlpfcRiboZero_Stats)), ],
					 CMC_Stats[match(ids,rownames(CMC_Stats)), ],
					 Caudate_Stats[match(ids,rownames(Caudate_Stats)), ],
					 Hippo_Stats[match(ids,rownames(Hippo_Stats)), ] )
rownames(mergedStats) = ids					 

sig_mergedStats = mergedStats[ rowSums( mergedStats[,grep("adj.P.Val",colnames(mergedStats))] < 0.05,na.rm=T )>0 ,	]				 

save(mergedStats, sig_mergedStats, file='/dcl01/lieber/ajaffe/Steve/SNX19/DE/rdas/differential_expression_results.rda')					 

write.csv(sig_mergedStats, '/dcl01/lieber/ajaffe/Steve/SNX19/DE/csvs/sig_mergedStats.csv',row.names=F)

colSums(sig_mergedStats[,grep("adj",colnames(sig_mergedStats),value=T)]	<0.05,na.rm=T)			 