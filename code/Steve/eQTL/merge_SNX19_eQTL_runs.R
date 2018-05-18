#qsub -l bluejay -l mf=30G,h_vmem=50G,h_stack=256M -M stephensemick@gmail.com -cwd -b y R CMD BATCH --no-save merge_SNX19_eQTL_runs.R

#Merge SNX19 eQTL runs
## Set up ##
library(jaffelab)
library(MatrixEQTL)
library(SummarizedExperiment)
library(matrixStats)
###########################
##  annotation object #####
###########################
#### DLPFC PolyA #### 
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_polyA/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_n495_updateMap_03082017.rda")
pd_DlpfcPolya = pd
#snp objects
snpMap_DlpfcPolya = snpMapAll
snpAll_DlpfcPolya = snpAll

#expression object
geneRpkm_DlpfcPolya = geneRpkm
geneMap_DlpfcPolya = geneMap

exonRpkm_DlpfcPolya = exonRpkm
exonMap_DlpfcPolya = exonMap

jRpkm_DlpfcPolya = jRpkm
jMap_DlpfcPolya = jMap

yExprs_DlpfcPolya = rbind(geneRpkm_DlpfcPolya, exonRpkm_DlpfcPolya, jRpkm_DlpfcPolya )
#### DLPFC RiboZero #### 
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_RiboZero/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_RiboZero_n485_updateMap.rda")
pd_DlpfcRibo = pd
#snp objects
snpMap_DlpfcRibo = snpMapAll
snpAll_DlpfcRibo = snpAll

#expression object
geneRpkm_DlpfcRibo= geneRpkm
geneMap_DlpfcRibo = geneMap

exonRpkm_DlpfcRibo = exonRpkm
exonMap_DlpfcRibo = exonMap

jRpkm_DlpfcRibo = jRpkm
jMap_DlpfcRibo = jMap

yExprs_DlpfcRibo = rbind(geneRpkm_DlpfcRibo, exonRpkm_DlpfcRibo, jRpkm_DlpfcRibo)
#### Caudate #### 
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/Caudate_RiboZero/rawAndRpkmCounts_plusGenotype_SNX19_Caudate_n468_updateMap.rda")
pd_Caudate = pd
#snp objects
snpMap_Caudate = snpMapAll
snpAll_Caudate = snpAll

#expression object
geneRpkm_Caudate = geneRpkm
geneMap_Caudate = geneMap

exonRpkm_Caudate = exonRpkm
exonMap_Caudate = exonMap

jRpkm_Caudate = jRpkm
jMap_Caudate = jMap

yExprs_Caudate = rbind(geneRpkm_Caudate, exonRpkm_Caudate, jRpkm_Caudate)
#### Hippo #### 
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/HIPPO_RiboZero/rawAndRpkmCounts_plusGenotype_SNX19_HIPPO_RiboZero_n452_updateMap.rda")
pd_Hippo = pd
#snp objects
snpMap_Hippo = snpMapAll
snpAll_Hippo= snpAll

#expression object
geneRpkm_Hippo = geneRpkm
geneMap_Hippo = geneMap

exonRpkm_Hippo = exonRpkm
exonMap_Hippo = exonMap

jRpkm_Hippo = jRpkm
jMap_Hippo = jMap

yExprs_Hippo = rbind(geneRpkm_Hippo, exonRpkm_Hippo, jRpkm_Hippo)

#### CMC #### 
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/CMC/rawAndRpkmCounts_plusGenotype_SNX19_CMC_n547.rda")
## filter samples to postnatal
pd$Age_of_Death[pd$Age_of_Death=="90+"] <- NA
pd <- plyr::rename(pd, c("Age_of_Death"="age"))
pd$age <- as.numeric(pd$age)
pd$Race = plyr::mapvalues(pd$Ethnicity, c("African-American", "Caucasian"), c("AA","CAUC") )
pd$Dx = plyr::mapvalues(pd$Dx, c("SCZ"), c("Schizo") )
pd_CMC= pd
#snp objects
snpMap_CMC = snpMapAll
snpAll_CMC = snpMatAll

#expression object
geneRpkm_CMC = geneRpkm
geneMap_CMC = geneMap

exonRpkm_CMC = exonRpkm
exonMap_CMC = exonMap

jRpkm_CMC = jRpkm
jMap_CMC = jMap

yExprs_CMC = rbind(geneRpkm_CMC, exonRpkm_CMC, jRpkm_CMC)

### Merge expression annotation object ###
geneMap = unique(rbind(geneMap_DlpfcPolya, geneMap_DlpfcRibo, geneMap_Caudate, geneMap_Hippo ))
exonMap = unique(rbind(exonMap_DlpfcPolya, exonMap_DlpfcRibo, exonMap_Caudate, exonMap_Hippo ))
jMap = unique(c(jMap_DlpfcPolya, jMap_DlpfcRibo, jMap_Caudate, jMap_Hippo ))
jMap$newName = names(jMap)
jMap = as.data.frame(jMap)
## common expression map to annotate
geneMap$EnsemblGeneID = rownames(geneMap)
geneMap$Class = "InEns"
geneMap$Type = "Gene"

colnames(exonMap)[1] = "EnsemblGeneID"
exonMap$Class = "InEns"
exonMap$Type = "Exon"

colnames(jMap)[c(16,19,20)] = c(
			"Class", "EnsemblGeneID", "Symbol")
jMap <- plyr::rename(jMap,  #renaming some of the variable names in the posJxn object
                       c("seqnames"="Chr", 
                         'start' = 'Start', 
                         'end'='End')) #rename "chrpos" to "snp_coord"			
jMap$Type = "Junction"

name = c("Chr","Start","End","EnsemblGeneID", "Symbol", "Type", "Class")
exprsMap = rbind(as.data.frame(geneMap)[,name],
	as.data.frame(exonMap)[,name],
	as.data.frame(jMap)[,name])
exprsMap = DataFrame(exprsMap)
#Convert hg19 eIDs to ENSEMBL exon IDs
ensemblGTF = rtracklayer::import(con="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ensembl/Homo_sapiens.GRCh37.75.gtf", format="gtf")
ensEXONS = as.data.frame(ensemblGTF)[which(ensemblGTF$type=="exon"),c("seqnames","start","end","strand","gene_id","exon_id")]
names(ensEXONS) = c("Chr","Start","End","Strand","ensemblID","exon_ensemblID")
rownames(ensEXONS) = NULL
ensEXONS$eID = paste0("e", rownames(ensEXONS))
ensEXONS$SNX19_ID = paste0("chr", ensEXONS$Chr,":",ensEXONS$Start,"-", ensEXONS$End, "(",ensEXONS$Strand, ")")
exprsMap$txMatchID = rownames(exprsMap)
eeMatch = match(exprsMap$txMatchID,ensEXONS$SNX19_ID)
exprsMap$txMatchID[!is.na(eeMatch) ] = ensEXONS[eeMatch[!is.na(eeMatch)], 'eID']
## add Tx data ##
xx=load("/dcl01/lieber/ajaffe/Steve/SNX19/Data/ensembl_V75_feature_to_Tx.rda")
tx = CharacterList(vector("list", nrow(exprsMap)))

## just exons ##
mmExon = match(exprsMap$txMatchID, names(allTx))
tx[!is.na(mmExon)] = allTx[mmExon[!is.na(mmExon)]]
tx[exprsMap$Type == "Junction"] = jMap$ensemblTx
exprsMap$NumTx = elementNROWS(tx)
exprsMap$WhichTx = sapply(tx, paste0, collapse=";")
exprsMap$WhichTx[exprsMap$Type == "Gene"] = NA
exprsMap$MergeID = c(rownames(geneMap), rownames(exonMap), jMap$newName)
exprsMap = as.data.frame(exprsMap)

##### Merge snp objects ######
col_keep = c('CHR','SNP','POS','COUNTED','ALT','chrpos','Type','newRef','newCount')
snpMapMerge = rbind(snpMap_DlpfcPolya[, col_keep],
					snpMap_DlpfcRibo[, col_keep], 
					snpMap_Caudate[, col_keep],
					snpMap_Hippo[, col_keep] )
snpMapMerge = snpMapMerge[!duplicated(snpMapMerge),]
snpMapMerge$name = snpMap_DlpfcPolya[match(snpMapMerge$SNP, snpMap_DlpfcPolya$SNP),'name']

## Load PGC Data ##
load('/dcl01/lieber/ajaffe/Steve/SNX19/PGC_GWAS/SCZ52_may13_9444230snps-scz2-snp-results-ckqny-scz2snpres-rall.rda')
PGC_52 = PGC_52[PGC_52$CHR=="chr11",]
PGC_52$chrpos = paste0(PGC_52$CHR, ":", PGC_52$BP)
PGC_52 = PGC_52[!(duplicated(PGC_52$chrpos) | duplicated(PGC_52$chrpos, fromLast=TRUE) ), ] #drop both SNPs from duplicated positions

#################
## load eqtls ###
#################

## DLPFC polyA ##
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/DLPFC_polyA/rdas/rawEqtls_DLPFC_polyA_age13_merged.rda')
cisDlpfcPolya = eqtl
colnames(cisDlpfcPolya)[-(1:3)] = paste0("DLPFC_PolyA_" ,colnames(cisDlpfcPolya)[-(1:3)] )

## DLPFC RiboZero ##
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/DLPFC_RiboZero/rdas/rawEqtls_DLPFC_RiboZero_age13_merged.rda')
cisDlpfcRibo = eqtl
colnames(cisDlpfcRibo)[-(1:3)] = paste0("DLPFC_RiboZero_" ,colnames(cisDlpfcRibo)[-(1:3)] )

## Caudate ##
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Caudate/rdas/rawEqtls_Caudate_RiboZero_age13_merged.rda')
cisCaudate = eqtl
colnames(cisCaudate)[-(1:3)] = paste0("Caudate_" ,colnames(cisCaudate)[-(1:3)] )

## Hippo ##
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Hippo/rdas/rawEqtls_Hippo_RiboZero_age13_merged.rda')
cisHippo = eqtl
colnames(cisHippo)[-(1:3)] = paste0("Hippo_" ,colnames(cisHippo)[-(1:3)] )

#################
## Merging ######
#################

eqtlNames = unique(c(rownames(cisDlpfcPolya), rownames(cisDlpfcRibo),rownames(cisCaudate), rownames(cisHippo) ))

cisMerge = cbind(cisDlpfcPolya[match(eqtlNames, rownames(cisDlpfcPolya)),-(1:3)],
	cisDlpfcRibo[match(eqtlNames, rownames(cisDlpfcRibo)),-(1:3)],
	cisCaudate[match(eqtlNames, rownames(cisCaudate)),-(1:3)],
	cisHippo[match(eqtlNames, rownames(cisHippo)),-(1:3)])
	
cisMerge$SNP = ss(eqtlNames, "\\.", 1)
cisMerge$Feature = ss(eqtlNames, "\\.", 2)
cisMerge$UniqueID = eqtlNames

rownames(cisMerge)=NULL
cisMerge = cisMerge[,c('SNP', 'Feature', 'UniqueID', colnames(cisMerge)[!colnames(cisMerge) %in% c('SNP','Feature', 'UniqueID' )]) ]

## which triallelic SNPs to drop and those that don't have GWAS
tt = table(snpMapMerge$SNP) #SNPs with same name
yy = table(snpMapMerge$chrpos) #SNPs at identical position
dropSnps = unique(c(names(tt[tt>1]), names(yy[yy>1]) ) )
cisMerge = cisMerge[!cisMerge$SNP%in%dropSnps,]

## Annotating SNPs

## add snp annotation ##
mmSnp = match(cisMerge$SNP, snpMapMerge$SNP)
cisMerge$snpName = snpMapMerge$name[mmSnp]
cisMerge$snpChr = paste0("chr",snpMapMerge$CHR[mmSnp])
cisMerge$snpPos = snpMapMerge$POS[mmSnp]
cisMerge$snpChrPos = snpMapMerge$chrpos[mmSnp]
cisMerge$snpType = snpMapMerge$Type[mmSnp]
cisMerge$snpRef =  snpMapMerge$newRef[mmSnp]
cisMerge$snpCounted =  snpMapMerge$newCount[mmSnp]

## add PGC info ##
mmPGC = match(cisMerge$snpChrPos, PGC_52$chrpos)
cisMerge$pgc52_snpName = PGC_52$SNP[mmPGC]
cisMerge$pgc52_A1 = PGC_52$A1[mmPGC]
cisMerge$pgc52_A2 = PGC_52$A2[mmPGC]
cisMerge$pgc52_OR = PGC_52$OR[mmPGC]
cisMerge$pgc52_SE = PGC_52$SE[mmPGC]
cisMerge$pgc52_P = PGC_52$P[mmPGC]

## adding risk directionality ##
cisMerge$pgc52_riskAllele = ifelse(cisMerge$pgc52_OR > 1, 	cisMerge$pgc52_A1, cisMerge$pgc52_A2)				

## Add CMC ##
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/CMC/rdas/rawEqtls_CMC_DLPFC_age13_merged.rda')
cisCMC = eqtl
cisCMC$snpChrPos = snpMap_CMC[match(cisCMC$SNP, snpMap_CMC$SNP),'chrpos']

mmCMC = match(paste0(cisMerge$snpChrPos,".",cisMerge$Feature), paste0(cisCMC$snpChrPos,".",cisCMC$Feature))
colnames(cisCMC)[-c(1:3,40)] = paste0("CMC_DLPFC_",colnames(cisCMC)[-c(1:3,40)]    )
cisMerge = cbind(cisMerge, cisCMC[mmCMC,-c(1:3,40)])
## Add exprs information ##
mmExprs = match(cisMerge$Feature, exprsMap$MergeID)
cisMerge = cbind(cisMerge, exprsMap[mmExprs,])
cisMerge$MergeID=NULL

#cisMerge$WhichTx = exprsMap[mmExprs,'WhichTx']
#cisMerge$NumTx = exprsMap[mmExprs,'NumTx']

#
## add ucsc ##
r = c(makeGRangesFromDataFrame(geneMap), makeGRangesFromDataFrame(exonMap), makeGRangesFromDataFrame(jMap))
exprsMap$BED = paste0(seqnames(r), ":", start(r), "-", end(r))
cisMerge$ucsc_feature = exprsMap$BED[match(cisMerge$Feature, exprsMap$MergeID)]

## snp ##
snpMapMerge$BED = paste0("chr", snpMapMerge$CHR, ":",
	snpMapMerge$POS -1,  "-", snpMapMerge$POS)
cisMerge$ucsc_snp = snpMapMerge$BED[match(cisMerge$SNP, snpMapMerge$SNP)]

colnames(cisMerge)[colnames(cisMerge) %in% c("Chr","Start","End","Type","Class")] = paste0("exprs",colnames(cisMerge)[colnames(cisMerge) %in% c("Chr","Start","End","Type","Class")])
save(cisMerge,  file='/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/Annotated_5LIBD-CMC_set_Merged.rda' )

## Filtering ##
cisMerge = cisMerge[!is.na(cisMerge$pgc52_P),] # Removing eQTLs where there is no PGC info for SNP
cisMerge = cisMerge[rowSums(is.na(cisMerge[,c("DLPFC_PolyA_Linear_All_statistic", "DLPFC_RiboZero_Linear_All_statistic", "Caudate_Linear_All_statistic", "Hippo_Linear_All_statistic")]) )==0,] #Removing eQTLs which are not tested in all LIBD datasets

## distance to feature ##
dToStart = cisMerge$snpPos - cisMerge$exprsStart #distance from single nucleotide polymorphism to expression start site
dToEnd = cisMerge$snpPos - cisMerge$exprsEnd #distance from single nucleotide polymorphism to expression end site
cisMerge$distSnpToFeature = apply(cbind(dToStart,dToEnd), 1, function(x) x[which.min(abs(x))]) #minimum distance

save(cisMerge,  file='/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing.rda' )
### Sig filter
sigFilter = rowSums(cisMerge[,grep("_pvalue", colnames(cisMerge))] < 1e-5, na.rm=TRUE)
cisMerge = cisMerge[sigFilter > 0,]

###### Add mean expression #######
## DLPFC PolyA 
yExprs_DlpfcPolya = yExprs_DlpfcPolya[match(cisMerge$Feature, rownames(yExprs_DlpfcPolya)), which(pd_DlpfcPolya$age>13 & pd_DlpfcPolya$Race %in% c("CAUC", "AA") ) ] 

## DLPFC RiboZero 
yExprs_DlpfcRibo = yExprs_DlpfcRibo[match(cisMerge$Feature, rownames(yExprs_DlpfcRibo)), which(pd_DlpfcRibo$age>13 & pd_DlpfcRibo$Race %in% c("CAUC", "AA") ) ] 

## Caudate
yExprs_Caudate = yExprs_Caudate[match(cisMerge$Feature, rownames(yExprs_Caudate)), which(pd_Caudate$age>13 & pd_Caudate$Race %in% c("CAUC", "AA") ) ] 

## Hippo
yExprs_Hippo = yExprs_Hippo[match(cisMerge$Feature, rownames(yExprs_Hippo)), which(pd_Hippo$Age>13 & pd_Hippo$Race %in% c("CAUC", "AA") ) ] 

## CMC
yExprs_CMC = yExprs_CMC[match(cisMerge$Feature, rownames(yExprs_CMC)), which(pd_CMC$age>13 & pd_CMC$Race %in% c("CAUC", "AA") ) ] 

##
yExprsList = list(DLPFC_polyA = yExprs_DlpfcPolya,
	DLPFC_RiboZero = yExprs_DlpfcRibo,
	Caudate_RiboZero = yExprs_Caudate,
	Hippo_RiboZero = yExprs_Hippo, 
	CMC_DLPFC = yExprs_CMC)
	
meanMat = sapply(yExprsList, function(x) rowMeans(x) )	
colnames(meanMat) = paste0("meanExprs_",colnames(meanMat))

cisMerge = cbind(cisMerge, meanMat)

###### Add SNP Frequency Information #######
## DLPFC PolyA 
xSNP_DlpfcPolya = snpAll_DlpfcPolya[match(cisMerge$SNP, rownames(snpAll_DlpfcPolya)), which(pd_DlpfcPolya$age>13 & pd_DlpfcPolya$Race %in% c("CAUC", "AA") ) ] 

## DLPFC RiboZero 
xSNP_DlpfcRibo = snpAll_DlpfcRibo[match(cisMerge$SNP, rownames(snpAll_DlpfcRibo)), which(pd_DlpfcRibo$age>13 & pd_DlpfcRibo$Race %in% c("CAUC", "AA") ) ] 

## Caudate
xSNP_Caudate = snpAll_Caudate[match(cisMerge$SNP, rownames(snpAll_Caudate)), which(pd_Caudate$age>13 & pd_Caudate$Race %in% c("CAUC", "AA") ) ] 

## Hippo
xSNP_Hippo = snpAll_Hippo[match(cisMerge$SNP, rownames(snpAll_Hippo)), which(pd_Hippo$Age>13 & pd_Hippo$Race %in% c("CAUC", "AA") ) ] 

## CMC
xSNP_CMC = snpAll_CMC[ snpMap_CMC[match(cisMerge$snpChrPos, snpMap_CMC$chrpos),'SNP'], which(pd_CMC$age>13 & pd_CMC$Race %in% c("CAUC", "AA") ) ] 

##
xSNPList = list(DLPFC_polyA = xSNP_DlpfcPolya,
	DLPFC_RiboZero = xSNP_DlpfcRibo,
	Caudate_RiboZero = xSNP_Caudate,
	Hippo_RiboZero = xSNP_Hippo,
	CMC_DLPFC = xSNP_CMC)
	
mafMat = sapply(xSNPList, function(x) rowSums(x, na.rm=TRUE)/(2*rowSums(!is.na(x))) )	
colnames(mafMat) = paste0("countedAlleleFreq_",colnames(mafMat))
cisMerge = cbind(cisMerge, mafMat)
########### Adding directionality #############

#### directionality
dataSets =c("DlpfcPolya","DlpfcRibo","Caudate", "Hippo", "CMC")
eqtlTraits =c("All",
			  "Control_AA",
			  "Control_CAUC",
			  "Schizo_AA",
			  "Schizo_CAUC",
			  "CAUC",
			  "AA")
length_eqtl_traits=length(eqtlTraits)
length_gwas_traits = length(dataSets)	
## set up matrics
snpMatchMat = eqtlSignMat = riskSignMat = matrix(0, 
	ncol = length_gwas_traits*length_eqtl_traits, nrow = nrow(cisMerge))
colnames(snpMatchMat) = paste0(rep(dataSets,
	each=length_eqtl_traits), "_", rep(eqtlTraits, times=length_gwas_traits))
colnames(eqtlSignMat) = colnames(riskSignMat) = colnames(snpMatchMat)

# does snp match?
snpMatchMat[cisMerge$pgc52_A2 == cisMerge$snpCounted,] = 1
snpMatchMat[grepl("I", cisMerge$pgc52_A2) & cisMerge$snpType == "Insertion",] = 1
snpMatchMat[cisMerge$pgc52_A2 == "D" & cisMerge$snpType == "Deletion",] = 1
snpMatchMat[cisMerge$pgc52_A1 == cisMerge$snpCounted,] = -1
snpMatchMat[grepl("I", cisMerge$pgc52_A2)& cisMerge$snpType == "Deletion",] = -1
snpMatchMat[cisMerge$pgc52_A2 == "D" & cisMerge$snpType == "Insertion",] = -1

# risk allele
riskSignMat[cisMerge$pgc52_A2 == cisMerge$pgc52_riskAllele, 1:35] = 1
riskSignMat[cisMerge$pgc52_A1 == cisMerge$pgc52_riskAllele, 1:35] = -1

# risk allele

#DLPFC PolyA
eqtlSignMat[, 1] =sign(cisMerge$DLPFC_PolyA_Linear_All_statistic)
eqtlSignMat[, 2] =sign(cisMerge$DLPFC_PolyA_Linear_Control_AA_statistic)
eqtlSignMat[, 3] =sign(cisMerge$DLPFC_PolyA_Linear_Control_CAUC_statistic)
eqtlSignMat[, 4] =sign(cisMerge$DLPFC_PolyA_Linear_Schizo_AA_statistic)
eqtlSignMat[, 5] =sign(cisMerge$DLPFC_PolyA_Linear_Schizo_CAUC_statistic)
eqtlSignMat[, 6] =sign(cisMerge$DLPFC_PolyA_Linear_CAUC_statistic)
eqtlSignMat[, 7] =sign(cisMerge$DLPFC_PolyA_Linear_AA_statistic)

#DLPFC RiboZero
eqtlSignMat[, 8] =sign(cisMerge$DLPFC_RiboZero_Linear_All_statistic)
eqtlSignMat[, 9] =sign(cisMerge$DLPFC_RiboZero_Linear_Control_AA_statistic)
eqtlSignMat[, 10] =sign(cisMerge$DLPFC_RiboZero_Linear_Control_CAUC_statistic)
eqtlSignMat[, 11] =sign(cisMerge$DLPFC_RiboZero_Linear_Schizo_AA_statistic)
eqtlSignMat[, 12] =sign(cisMerge$DLPFC_RiboZero_Linear_Schizo_CAUC_statistic)
eqtlSignMat[, 13] =sign(cisMerge$DLPFC_RiboZero_Linear_CAUC_statistic)
eqtlSignMat[, 14] =sign(cisMerge$DLPFC_RiboZero_Linear_AA_statistic)

# Caudate
eqtlSignMat[, 15] =sign(cisMerge$Caudate_Linear_All_statistic)
eqtlSignMat[, 16] =sign(cisMerge$Caudate_Linear_Control_AA_statistic)
eqtlSignMat[, 17] =sign(cisMerge$Caudate_Linear_Control_CAUC_statistic)
eqtlSignMat[, 18] =sign(cisMerge$Caudate_Linear_Schizo_AA_statistic)
eqtlSignMat[, 19] =sign(cisMerge$Caudate_Linear_Schizo_CAUC_statistic)
eqtlSignMat[, 20] =sign(cisMerge$Caudate_Linear_CAUC_statistic)
eqtlSignMat[, 21] =sign(cisMerge$Caudate_Linear_AA_statistic)

# Hippo
eqtlSignMat[, 22] =sign(cisMerge$Hippo_Linear_All_statistic)
eqtlSignMat[, 23] =sign(cisMerge$Hippo_Linear_Control_AA_statistic)
eqtlSignMat[, 24] =sign(cisMerge$Hippo_Linear_Control_CAUC_statistic)
eqtlSignMat[, 25] =sign(cisMerge$Hippo_Linear_Schizo_AA_statistic)
eqtlSignMat[, 26] =sign(cisMerge$Hippo_Linear_Schizo_CAUC_statistic)
eqtlSignMat[, 27] =sign(cisMerge$Hippo_Linear_CAUC_statistic)
eqtlSignMat[, 28] =sign(cisMerge$Hippo_Linear_AA_statistic)

# CMC DLPFC
eqtlSignMat[, 29] =sign(cisMerge$CMC_DLPFC_Linear_All_statistic)
eqtlSignMat[, 30] =sign(cisMerge$CMC_DLPFC_Linear_Control_AA_statistic)
eqtlSignMat[, 31] =sign(cisMerge$CMC_DLPFC_Linear_Control_CAUC_statistic)
eqtlSignMat[, 32] =sign(cisMerge$CMC_DLPFC_Linear_Schizo_AA_statistic)
eqtlSignMat[, 33] =sign(cisMerge$CMC_DLPFC_Linear_Schizo_CAUC_statistic)
eqtlSignMat[, 34] =sign(cisMerge$CMC_DLPFC_Linear_CAUC_statistic)
eqtlSignMat[, 35] =sign(cisMerge$CMC_DLPFC_Linear_AA_statistic)

# combine	
directionalityMat = snpMatchMat*riskSignMat*eqtlSignMat
directionalityMat[directionalityMat==1] = "U"
directionalityMat[directionalityMat==-1] = "D"
directionalityMat[is.na(directionalityMat)] = "M"

# make strings
dirString = apply(directionalityMat,1, function(x) {
	paste(paste(x[1:7], collapse=""), 
		paste(x[8:14], collapse=""),  
		paste(x[15:21], collapse=""),
		paste(x[22:28], collapse=""),
		paste(x[29:35], collapse=""),
		sep="_")
})
cisMerge$riskDirectString = dirString

########## Saving final object ###############
save(cisMerge,  file='/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt.rda')
write.csv(cisMerge,row.names=FALSE, file = '/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/tables/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt.csv')
write.csv(cisMerge,row.names=FALSE, file = gzfile('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/tables/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt.csv.gz') )

########## Filter down based on linear results only ##############
### Sig filter
#sigFilter = rowSums(cisMerge[,grepl("_pvalue", colnames(cisMerge)) & grepl("Linear", colnames(cisMerge))] < 1e-4, na.rm=TRUE)
#cisMerge = cisMerge[sigFilter > 0,]


### Additional annotation
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt.rda')

## Load PGC 2014 ##
load('/dcl01/lieber/ajaffe/Steve/SNX19/PGC_GWAS/PGC_2014.rda')
PGC_2014 = PGC_2014[PGC_2014$CHR=="11",]
PGC_2014$CHR = paste0("chr",PGC_2014$CHR )
PGC_2014$chrpos = paste0(PGC_2014$CHR, ":", PGC_2014$BP)
PGC_2014 = PGC_2014[!(duplicated(PGC_2014$chrpos) | duplicated(PGC_2014$chrpos, fromLast=TRUE) ), ] #drop both SNPs from duplicated positions

## add PGC info 2014 ##
mmPGC2014 = match(cisMerge$snpChrPos, PGC_2014$chrpos)
cisMerge$pgc2014_snpName = PGC_2014$SNP[mmPGC2014]
cisMerge$pgc2014_A1 = PGC_2014$A1[mmPGC2014]
cisMerge$pgc2014_A2 = PGC_2014$A2[mmPGC2014]
cisMerge$pgc2014_OR = PGC_2014$OR[mmPGC2014]
cisMerge$pgc2014_SE = PGC_2014$SE[mmPGC2014]
cisMerge$pgc2014_P = PGC_2014$P[mmPGC2014]

## Load PGC 2016 ##
load('/dcl01/lieber/ajaffe/Steve/SNX19/PGC_GWAS/PGC2016_daner_PGC_SCZ_1016_102716_8736907snps.rda')
PGC_2016 = PGC_2016[PGC_2016$CHR=="11",]
PGC_2016$CHR = paste0("chr",PGC_2016$CHR )
PGC_2016$chrpos = paste0(PGC_2016$CHR, ":", PGC_2016$BP)
PGC_2016 = PGC_2016[!(duplicated(PGC_2016$chrpos) | duplicated(PGC_2016$chrpos, fromLast=TRUE) ), ] #drop both SNPs from duplicated positions

## add PGC info 2016 ##
mmPGC2016 = match(cisMerge$snpChrPos, PGC_2016$chrpos)
cisMerge$pgc2016_snpName = PGC_2016$SNP[mmPGC2016]
cisMerge$pgc2016_A1 = PGC_2016$A1[mmPGC2016]
cisMerge$pgc2016_A2 = PGC_2016$A2[mmPGC2016]
cisMerge$pgc2016_OR = PGC_2016$OR[mmPGC2016]
cisMerge$pgc2016_SE = PGC_2016$SE[mmPGC2016]
cisMerge$pgc2016_P = PGC_2016$P[mmPGC2016]

#load and move fields around
colOrder = readxl::read_excel('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/SNX19_FieldNames-cur2.xlsx')
colOrder$FIELDNAME[!colOrder$FIELDNAME%in%colnames(cisMerge)]
cisMerge = cisMerge[,colOrder$FIELDNAME]
library(jaffelab)
cisMerge$riskDirectString = paste(ss(cisMerge$riskDirectString,"_" ,1),  ss(cisMerge$riskDirectString,"_" ,2), ss(cisMerge$riskDirectString,"_" ,5), ss(cisMerge$riskDirectString,"_" ,4), ss(cisMerge$riskDirectString,"_" ,3), sep="_" )
########## Saving final object ###############
save(cisMerge,  file='/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC.rda')
write.csv(cisMerge,row.names=FALSE, file = '/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/tables/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC.csv')
write.csv(cisMerge,row.names=FALSE, file = gzfile('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/tables/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC.csv.gz') )
