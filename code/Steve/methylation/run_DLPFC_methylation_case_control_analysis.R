###
library(jaffelab)
library(MatrixEQTL)
library(SummarizedExperiment)
library(matrixStats)
setwd('/dcl01/lieber/ajaffe/Steve/SNX19/')
## load counts
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_polyA/SNX19_DNAm_data_DLPFC.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_polyA/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_n495_updateMap_03082017.rda")

## reducing to shared people and matching
kk = match(methPheno$BrNum, pd$BrNum)
methPheno = methPheno[!is.na(kk),]
methMat = methMat[,!is.na(kk)]
pd = pd[kk[!is.na(kk)], ]
snpAll = snpAll[,pd$BrNum]

## filter samples to postnatal
keepIndex = which(pd$Race %in% c("CAUC", "AA") &  pd$age > 13)
methPheno = methPheno[keepIndex,]
methMat = methMat[,keepIndex]
pd = pd[keepIndex,]
snpAll = snpAll[,keepIndex]

#####
load("/dcl01/lieber/ajaffe/Steve/SNX19/meqtl_runs/DLPFC/rdas/meth_PCs_DLPFC_polyA_age13_matchedUp.rda")
mod = model.matrix(~methPheno$Dx + methPheno$Age + methPheno$Gender + methPheno$snpPC1 + methPheno$snpPC2 + methPheno$snpPC3 + methPheno$snpPC4 + methPheno$snpPC5 + pca_meth[methPheno$BrNum,1:15] )
fit = lmFit(methMat, mod)
fit <- eBayes(fit)
methDLPFC_stats = topTable(fit, coef='methPheno$DxSchizo', n=Inf, sort.by = 'p',genelist = )
write.csv(methDLPFC_stats, file='/dcl01/lieber/ajaffe/Steve/SNX19/tables/caseControl_differential_methylation_N401_442cpg.csv',row.names=T)
###########################
### get expression PCs ####
#pca_meth = prcomp(t(methMat)) #principle component decomposition on the transpose of that matrix + 1
#pca_meth = pca_meth$x[,1:15] #take the first 15 PCs from that matrix
#
#########################
### set up model matrix
#pd$Dx = factor(pd$Dx, levels=c("Control","Schizo") )
#mod = model.matrix(~pd$snpPC1 + pd$snpPC2 + pd$snpPC3 + pd$snpPC4 + pd$snpPC5 + pca_meth) #forming the model matrix from the SNP PCs and the expression PCs
#colnames(mod)[2:ncol(mod)] = c(paste0("snpPC",1:5), paste0("methPC",1:15)) #Renaming the column names of the model matrix
#covs = SlicedData$new(t(mod[,-1])) #This part employs the "MatrixEQTL" package
#
#### set up SNPs
#snpSlice = SlicedData$new(as.matrix(snpAll) ) #formating the snp data for the MatrixEQTL package 
#snpSlice$ResliceCombined(sliceSize = 5000) 
#snpspos = snpMapAll[,c("SNP","CHR","POS")]
#snpspos$CHR = paste0("chr",snpspos$CHR) #concatenating the string "chr" with the numbers in "CHR"
#colnames(snpspos) = c("SNP","chr","pos")
#rownames(snpspos) = NULL
#
#### set up expression data
#yMeth = methMat
#meth = SlicedData$new(as.matrix(yMeth))
#meth$ResliceCombined(sliceSize = 5000)
#
### make position data frame
#posMeth = methMap[,c("chr","pos","pos")]
#colnames(posMeth) = c("Chr","Start","End")
#posMeth$ProbeID = rownames(methMap)
#posMeth = posMeth[,c(4,1:3)]
#
##### full eqtl run ####
# methQtl_all = Matrix_eQTL_main(snps=snpSlice, 
#	 gene = meth, 
#	 cvrt = covs, 
#	 output_file_name.cis =  "cis.txt" ,
#	 output_file_name =  "trans.txt" ,
#	 pvOutputThreshold.cis = 0,  
#	 pvOutputThreshold=1,
#	 snpspos = snpspos, 
#	 genepos = posMeth, 
#	 useModel = modelLINEAR,	
#	 cisDist=1e6,
#	 pvalue.hist = 100,
#	 min.pv.by.genesnp = TRUE)
#
##### eqtl ####
#meqtl = methQtl_all$all$eqtls #extracting cis eqtl information from the eQTL analysis
#colnames(meqtl)[1:2] = c("SNP","ProbeID") #adding a column called feature
#meqtl$ProbeID = as.character(meqtl$ProbeID) #the ProbeID column contains expression information
#meqtl$SNP = as.character(meqtl$SNP) #converting the eqtl snps variable to a character variable
#colnames(meqtl)[colnames(meqtl) %in% c('statistic','pvalue','FDR','beta')] = paste0("All_", colnames(meqtl)[colnames(meqtl) %in% c('statistic','pvalue','FDR','beta')]) #Adding "ALL_ to "statistics, p-value, FDR, and beta" to help distinguish from future analyses 
#meqtl$UniqueID = paste0(meqtl$SNP, ".", meqtl$ProbeID) #making a variable called Unique ID... Just a combination of SNP and feature
#meqtl = meqtl[,c('UniqueID', 'SNP', 'ProbeID', 'All_statistic', 'All_pvalue', 'All_FDR', 'All_beta')] #just reordering the columns of the dataframe	 
#	 
############################################
########### subgroups ###################
############################################
#
##### Splitting Into Subgroups ####
##Now we may be interested in eQTL analysis among a particular subgroup (e.g. African Americans)
#splitit <- function(x) split(seq(along=x),x) #custom function written to split categories into different subgroups
#
##for meQTL analysis across different subgroups (race and disease status)
#gIndexes = splitit(paste0(pd$Dx, "_", pd$Race)) #meqtl analysis for four subgroups: Caucasian controls, African American controls, Caucasian Schizophrenics, and African American Schizophrenics.
#gIndexes$CAUC = which(pd$Race == "CAUC") #meqtl analysis for Caucasians
#gIndexes$AA = which(pd$Race == "AA") #meqtl analysis for African Americans
#
##### meQTL Analysis for Each Subgroup Defined Above
##Now we just apply a loop to perform meqtl analysis for each group defined above. 
#gStatList = parallel::mclapply(gIndexes, function(ii) {
#  yy = SlicedData$new( as.matrix(yMeth[,ii]) )
#  ssnp = SlicedData$new( as.matrix(snpAll[,ii]) )
#  covs = SlicedData$new(t(mod[ii,-1]))
#  me = Matrix_eQTL_main(snps=ssnp, 
#                        gene = yy, 
#                        cvrt = covs,
#                        output_file_name.cis =  "cis.txt" ,
#						output_file_name =  "trans.txt" ,
#                        pvOutputThreshold.cis = 0,  
#                        pvOutputThreshold=1,
#                        snpspos = snpspos, 
#                        genepos = posMeth, 
#                        useModel = modelLINEAR,
#						cisDist=1e6,
#                        pvalue.hist = 100,
#                        min.pv.by.genesnp = TRUE)
#  e = me$all$eqtls
#  e$UniqueID = paste0(e$snps, ".", e$gene)
#  e = e[match(meqtl$UniqueID, e$UniqueID),]
#  rownames(e) = e$UniqueID
#  e[,c(3,4,6)]
#}, mc.cores=1) #mc cores corresponds to number of cores available for multicore processing
#gStats = do.call("cbind", gStatList) #combining all these different subgroup analyses into one
#colnames(gStats) = gsub("\\.", "_", colnames(gStats)) #replacing periods with underscores
#meqtl = cbind(meqtl, gStats) #add the subgroup meqtl analysis
#colnames(meqtl)[-(1:3)] = paste0("Linear_" , colnames(meqtl)[-(1:3)])
#
##### full meqtl run ANOVA ####
# methQtl_all_anova = Matrix_eQTL_main(snps=snpSlice, 
#	 gene = meth, 
#	 cvrt = covs, 
#	 output_file_name.cis =  "cis.txt" ,
#	 output_file_name =  "trans.txt" ,
#	 pvOutputThreshold.cis = 0,
#	 pvOutputThreshold=1,
#	 snpspos = snpspos,
#	 genepos = posMeth,
#	 useModel = modelANOVA,
#	 cisDist=1e6,
#	 pvalue.hist = 100,
#	 min.pv.by.genesnp = TRUE)
#	 
#meqtl_anova = methQtl_all_anova$all$eqtls
#meqtl_anova$UniqueID = paste0(meqtl_anova$snps, ".", meqtl_anova$gene)
#meqtl_anova = meqtl_anova[match(meqtl$UniqueID, meqtl_anova$UniqueID),]
#rownames(meqtl_anova) = meqtl_anova$UniqueID
#meqtl_anova = meqtl_anova[,3:4]
#colnames(meqtl_anova) = paste0( "All_", colnames(meqtl_anova)   )
############################################
########### anova subgroups ###################
############################################
#
##### Splitting Into Subgroups ####
#
##### eQTL Analysis for Each Subgroup Defined Above
##Now we just apply a loop to perform eQTL analysis for each group defined above. 
#gStatList = parallel::mclapply(gIndexes, function(ii) {
#  yy = SlicedData$new( as.matrix(yMeth[,ii]) )
#  ssnp = SlicedData$new( as.matrix(snpAll[,ii]) )
#  covs = SlicedData$new(t(mod[ii,-1]))
#  me = Matrix_eQTL_main(snps=ssnp, 
#                        gene = yy, 
#                        cvrt = covs,
#                        output_file_name.cis =  "cis.txt" ,
#						output_file_name =  "trans.txt" ,
#                        pvOutputThreshold.cis = 0,  
#                        pvOutputThreshold=1,
#                        snpspos = snpspos, 
#                        genepos = posMeth, 
#                        useModel = modelANOVA,
#						cisDist=1e6,
#                        pvalue.hist = 100,
#                        min.pv.by.genesnp = TRUE)
#  e = me$all$eqtls
#  e$UniqueID = paste0(e$snps, ".", e$gene)
#  e = e[match(meqtl$UniqueID, e$UniqueID),]
#  rownames(e) = e$UniqueID
#  e[,c(3,4,6)]
#}, mc.cores=1) #mc cores corresponds to number of cores available for multicore processing
#gStats = do.call("cbind", gStatList) #combining all these different subgroup analyses into one
#colnames(gStats) = gsub("\\.", "_", colnames(gStats)) #replacing periods with underscores
#
####
#meqtl_anova = cbind(meqtl_anova,gStats)
#colnames(meqtl_anova)[] = paste0("ANOVA_" , colnames(meqtl_anova)[])
#
#### Combining ANOVA and Linear model results and saving
#meqtl = cbind(meqtl, meqtl_anova) #add the subgroup eQTL analysis
#save(meqtl, file = '/dcl01/lieber/ajaffe/Steve/SNX19/meqtl_runs/DLPFC/rdas/rawMeqtls_DLPFC_polyA_age13_merged.rda')

### Annotate meQTL Results ### 
#Based on: /users/ajaffe/Lieber/Projects/450k/ECD2014/check_meQTLs_adult.R
load('/dcl01/lieber/ajaffe/Steve/SNX19/meqtl_runs/DLPFC/rdas/rawMeqtls_DLPFC_polyA_age13_merged.rda')
#load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_polyA/SNX19_DNAm_data_DLPFC.rda")
#load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_polyA/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_n495_updateMap_03082017.rda")
##

## annotate
m = match(meqtl$SNP, snpMapAll$SNP)
meqtl$snpChr = paste0("chr", snpMapAll$CHR[m])
meqtl$snpPos = snpMapAll$POS[m]
meqtl$snpChrPos = snpMapAll$chrpos[m]
meqtl$snpRsNum = snpMapAll$name[m]
meqtl$snpCounted = snpMapAll$newCount[m]
meqtl$snpAlt = snpMapAll$newRef[m]
meqtl$inSampleMAF = rowSums(snpAll[m,],na.rm=TRUE)/
	(2*rowSums(!is.na(snpAll[m,])))
	
### see if SNP distrupts a CpG
library(BSgenome.Hsapiens.UCSC.hg19)
gr = GRanges(meqtl$snpChr, IRanges(meqtl$snpPos-1, meqtl$snpPos+1))
trio = getSeq(Hsapiens, gr)

meqtl$disruptCpG = vcountPattern("CG", trio)
meqtl$methChr = as.character(methMap[match(meqtl$ProbeID, rownames(methMap)),'chr' ])
meqtl$methPos = as.numeric(methMap[match(meqtl$ProbeID, rownames(methMap)),'pos' ])
meqtl$methProbeRs = methMap[match(meqtl$ProbeID, rownames(methMap)),'Probe_rs' ]
meqtl$methCpgRs = methMap[match(meqtl$ProbeID, rownames(methMap)),'CpG_rs' ]

meqtl$distMethToSnp = meqtl$methPos - meqtl$snpPos
meqtl$meanMeth = rowMeans(methMat)[meqtl$ProbeID]

############# PGC Information ##################
## Load PGC Data ##
load('/dcl01/lieber/ajaffe/Steve/SNX19/PGC_GWAS/SCZ52_may13_9444230snps-scz2-snp-results-ckqny-scz2snpres-rall.rda')
PGC_52 = PGC_52[PGC_52$CHR=="chr11",]
PGC_52$chrpos = paste0(PGC_52$CHR, ":", PGC_52$BP)
PGC_52 = PGC_52[!(duplicated(PGC_52$chrpos) | duplicated(PGC_52$chrpos, fromLast=TRUE) ), ] #drop both SNPs from duplicated positions

## add 52 PGC info ##
mmPGC = match(meqtl$snpChrPos, PGC_52$chrpos)
meqtl$pgc52_snpName = PGC_52$SNP[mmPGC]
meqtl$pgc52_A1 = PGC_52$A1[mmPGC]
meqtl$pgc52_A2 = PGC_52$A2[mmPGC]
meqtl$pgc52_OR = PGC_52$OR[mmPGC]
meqtl$pgc52_SE = PGC_52$SE[mmPGC]
meqtl$pgc52_P = PGC_52$P[mmPGC]

## Load PGC 2014 ##
load('/dcl01/lieber/ajaffe/Steve/SNX19/PGC_GWAS/PGC_2014.rda')
PGC_2014 = PGC_2014[PGC_2014$CHR=="11",]
PGC_2014$CHR = paste0("chr",PGC_2014$CHR )
PGC_2014$chrpos = paste0(PGC_2014$CHR, ":", PGC_2014$BP)
PGC_2014 = PGC_2014[!(duplicated(PGC_2014$chrpos) | duplicated(PGC_2014$chrpos, fromLast=TRUE) ), ] #drop both SNPs from duplicated positions

## add PGC info 2014 ##
mmPGC2014 = match(meqtl$snpChrPos, PGC_2014$chrpos)
meqtl$pgc2014_snpName = PGC_2014$SNP[mmPGC2014]
meqtl$pgc2014_A1 = PGC_2014$A1[mmPGC2014]
meqtl$pgc2014_A2 = PGC_2014$A2[mmPGC2014]
meqtl$pgc2014_OR = PGC_2014$OR[mmPGC2014]
meqtl$pgc2014_SE = PGC_2014$SE[mmPGC2014]
meqtl$pgc2014_P = PGC_2014$P[mmPGC2014]

## Load PGC 2016 ##
load('/dcl01/lieber/ajaffe/Steve/SNX19/PGC_GWAS/PGC2016_daner_PGC_SCZ_1016_102716_8736907snps.rda')
PGC_2016 = PGC_2016[PGC_2016$CHR=="11",]
PGC_2016$CHR = paste0("chr",PGC_2016$CHR )
PGC_2016$chrpos = paste0(PGC_2016$CHR, ":", PGC_2016$BP)
PGC_2016 = PGC_2016[!(duplicated(PGC_2016$chrpos) | duplicated(PGC_2016$chrpos, fromLast=TRUE) ), ] #drop both SNPs from duplicated positions

## add PGC info 2016 ##
mmPGC2016 = match(meqtl$snpChrPos, PGC_2016$chrpos)
meqtl$pgc2016_snpName = PGC_2016$SNP[mmPGC2016]
meqtl$pgc2016_A1 = PGC_2016$A1[mmPGC2016]
meqtl$pgc2016_A2 = PGC_2016$A2[mmPGC2016]
meqtl$pgc2016_OR = PGC_2016$OR[mmPGC2016]
meqtl$pgc2016_SE = PGC_2016$SE[mmPGC2016]
meqtl$pgc2016_P = PGC_2016$P[mmPGC2016]

save(meqtl, file = '/dcl01/lieber/ajaffe/Steve/SNX19/meqtl_runs/DLPFC/rdas/annotatedMeqtls_DLPFC_polyA_age13_merged.rda')
#write.csv(meqtl, file = '/dcl01/lieber/ajaffe/Steve/SNX19/meqtl_runs/DLPFC/tables/annotatedMeqtls_DLPFC_polyA_age13_merged.csv')
#write.csv(meqtl, file = gzfile('/dcl01/lieber/ajaffe/Steve/SNX19/meqtl_runs/DLPFC/tables/annotatedMeqtls_DLPFC_polyA_age13_merged.csv.gz') )

### Cut after filtering to nominal significance
sigFilter = rowSums(meqtl[,grep("_pvalue", colnames(meqtl))] < 1e-3, na.rm=TRUE)
meqtl_cut = meqtl[sigFilter > 0,]

save(meqtl_cut, file = '/dcl01/lieber/ajaffe/Steve/SNX19/meqtl_runs/DLPFC/rdas/annotatedMeqtls_DLPFC_polyA_age13_merged_1e-3sigFilt.rda')
write.csv(meqtl_cut, file = '/dcl01/lieber/ajaffe/Steve/SNX19/meqtl_runs/DLPFC/tables/annotatedMeqtls_DLPFC_polyA_age13_merged_1e-3sigFilt.csv')
write.csv(meqtl_cut, file = gzfile('/dcl01/lieber/ajaffe/Steve/SNX19/meqtl_runs/DLPFC/tables/annotatedMeqtls_DLPFC_polyA_age13_merged_1e-3sigFilt.csv.gz') )


### add gene info
#library(bumphunter)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#theTranscripts = annotateTranscripts( TxDb.Hsapiens.UCSC.hg19.knownGene,codingOnly=TRUE)
#
#an = annotateNearest(map, theTranscripts)	
#map$nearestGene = as.character(theTranscripts$Gene)[an$subjectHits]
#map$nearestGeneDist = an$dist
#	
#meqtl$nearestGene = map$nearestGene[match(meqtl$cpg, names(map))]	
#meqtl$nearestGeneDist = map$nearestGeneDist[match(meqtl$cpg, names(map))]	
