#qsub -l bluejay -l mf=5G,h_vmem=15G,h_stack=256M -M stephensemick@gmail.com -cwd -b y R CMD BATCH --no-save run_Caudate_RiboZero_eqtls.R

###
library(jaffelab)
library(MatrixEQTL)
library(SummarizedExperiment)
library(matrixStats)
setwd('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Caudate')
## load counts
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/Caudate_RiboZero/rawAndRpkmCounts_plusGenotype_SNX19_Caudate_n468_updateMap.rda")

## filter samples to postnatal
keepIndex = which(pd$Race %in% c("CAUC", "AA") &  pd$age > 13)
pd = pd[keepIndex,]
geneRpkm = geneRpkm[,keepIndex]
exonRpkm = exonRpkm[,keepIndex]
jRpkm = jRpkm[,keepIndex]

### filter to people in common ###
mm = match(pd$BrNum, colnames(snpAll ))
geneRpkm = geneRpkm[,!is.na(mm)]
exonRpkm = exonRpkm[,!is.na(mm)]
jRpkm = jRpkm[,!is.na(mm)]

snpAll  = snpAll [,mm[!is.na(mm)]]

### filter features ###
gIndex=which(rowMeans(geneRpkm) > 0.01) #selecting out genes with sufficient expression (genes with less than 0.01 read per kilobase million are not perserved)
geneRpkm = geneRpkm[gIndex,] #keeping only the genes above the rpkm threshold (threshold 0.01 here)
geneMap = geneMap[gIndex,] #keeping only the gene map for genes above the rpkm threshold

eIndex=which(rowMeans(exonRpkm) > 0.01)  #selecting out exons with sufficient expression (exons with less than 0.01 read per kilobase million are not perserved)
exonRpkm = exonRpkm[eIndex,] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
exonMap = exonMap[eIndex,] #keeping only the gene map for genes above the rpkm threshold

jIndex = which(rowMeans(jRpkm) > 0.2)#There are an excessive number of junctions. Apply a stricter filter and drop novel junctions to reduce the number.
jRpkm = jRpkm[jIndex,] #keeping only the junctions above the rpkm threshold (threshold = 0.01 rpkm here) and are NOT novel
jMap= jMap[jIndex] #keeping only the junction map for junctions above the rpkm threshold and are NOT novel

##########################
## get expression PCs ####
exprs3features = log2(rbind(geneRpkm, exonRpkm, jRpkm)+1) #combine the three types of data into a single matrix
pcaexprs3features = prcomp(t(exprs3features+1)) #principle component decomposition on the transpose of that matrix + 1
all3PCs = pcaexprs3features$x[,1:15] #take the first 15 PCs from that matrix
save(all3PCs, file="rdas/all3Feature_PCs_Caudate_RiboZero_age13_matchedUp.rda")

########################
## set up model matrix
pd$Dx = factor(pd$Dx, levels=c("Control","Schizo", "Bipolar") )
mod = model.matrix(~pd$snpPC1 + pd$snpPC2 + pd$snpPC3 + pd$snpPC4 + pd$snpPC5 + all3PCs) #forming the model matrix from the SNP PCs and the expression PCs
colnames(mod)[2:ncol(mod)] = c(paste0("snpPC",1:5), paste0("exprsPC",1:15)) #Renaming the column names of the model matrix
covs = SlicedData$new(t(mod[,-1])) #This part employs the "MatrixEQTL" package

### set up SNPs
snpSlice = SlicedData$new(as.matrix(snpAll) ) #formating the snp data for the MatrixEQTL package 
snpSlice$ResliceCombined(sliceSize = 5000) 
snpspos = snpMapAll[,c("SNP","CHR","POS")]
snpspos$CHR = paste0("chr",snpspos$CHR) #concatenating the string "chr" with the numbers in "CHR"
colnames(snpspos) = c("SNP","chr","pos")
rownames(snpspos) = NULL

### set up expression data
yExprs = exprs3features
exprs = SlicedData$new(as.matrix(yExprs))
exprs$ResliceCombined(sliceSize = 5000)

## make position data frame
posGene = geneMap[,c('Chr', 'Start', 'End')] #gene position information taken from the geneMap
posGene$name = rownames(geneMap) #gene position names variable made from rownames of geneMap
posExon = exonMap[,c('Chr', 'Start', 'End')] #exon position information taken from exonMap
posExon$name = rownames(exonMap) #exon position names variable made from rownames of exonMap
posJxn = as.data.frame(jMap)[,c('seqnames','start','end')] #extracting junction position information from the jMap object after converting it to a dataframe
colnames(posJxn) = c("Chr","Start","End")
posJxn$name = names(jMap)
posExprs = rbind(posGene, posExon, posJxn) #combining position information for the three different types of expression data
posExprs = posExprs[,c(4,1:3)]

#### full eqtl run ####
 meQtl_all = Matrix_eQTL_main(snps=snpSlice, 
	 gene = exprs, 
	 cvrt = covs, 
	 output_file_name.cis =  "cis.txt" ,
	 output_file_name =  "trans.txt" ,
	 pvOutputThreshold.cis = 0,  
	 pvOutputThreshold=1,
	 snpspos = snpspos, 
	 genepos = posExprs, 
	 useModel = modelLINEAR,	
	 cisDist=1e6,
	 pvalue.hist = 100,
	 min.pv.by.genesnp = TRUE)

#### eqtl ####
eqtl = meQtl_all$all$eqtls #extracting cis eqtl information from the eQTL analysis
colnames(eqtl)[1:2] = c("SNP","Feature") #adding a column called feature
eqtl$Feature = as.character(eqtl$Feature) #the feature column contains expression information
eqtl$SNP = as.character(eqtl$SNP) #converting the eqtl snps variable to a character variable
colnames(eqtl)[colnames(eqtl) %in% c('statistic','pvalue','FDR','beta')] = paste0("All_", colnames(eqtl)[colnames(eqtl) %in% c('statistic','pvalue','FDR','beta')]) #Adding "ALL_ to "statistics, p-value, FDR, and beta" to help distinguish from future analyses 
eqtl$UniqueID = paste0(eqtl$SNP, ".", eqtl$Feature) #making a variable called Unique ID... Just a combination of SNP and feature
eqtl = eqtl[,c('UniqueID', 'SNP', 'Feature', 'All_statistic', 'All_pvalue', 'All_FDR', 'All_beta')] #just reordering the columns of the dataframe	 
	 
###########################################
########## subgroups ###################
###########################################

#### Splitting Into Subgroups ####
#Now we may be interested in eQTL analysis among a particular subgroup (e.g. African Americans)
splitit <- function(x) split(seq(along=x),x) #custom function written to split categories into different subgroups

#for eQTL analysis across different subgroups (race and disease status)
gIndexes = splitit(paste0(pd$Dx, "_", pd$Race)) #eQTL analysis for four subgroups: Caucasian controls, African American controls, Caucasian Schizophrenics, and African American Schizophrenics.
gIndexes$CAUC = which(pd$Race == "CAUC") #eQTL analysis for Caucasians
gIndexes$AA = which(pd$Race == "AA") #eQTL analysis for African Americans
#DROPPING Bipolar_AA b/c n=4, sample size is too small for analysis
gIndexes<- gIndexes[!names(gIndexes)%in%"Bipolar_AA"]

#### eQTL Analysis for Each Subgroup Defined Above
#Now we just apply a loop to perform eQTL analysis for each group defined above. 
gStatList = parallel::mclapply(gIndexes, function(ii) {
  yy = SlicedData$new( as.matrix(yExprs[,ii]) )
  ssnp = SlicedData$new( as.matrix(snpAll[,ii]) )
  covs = SlicedData$new(t(mod[ii,-1]))
  me = Matrix_eQTL_main(snps=ssnp, 
                        gene = yy, 
                        cvrt = covs,
                        output_file_name.cis =  "cis.txt" ,
						output_file_name =  "trans.txt" ,
                        pvOutputThreshold.cis = 0,  
                        pvOutputThreshold=1,
                        snpspos = snpspos, 
                        genepos = posExprs, 
                        useModel = modelLINEAR,
						cisDist=1e6,
                        pvalue.hist = 100,
                        min.pv.by.genesnp = TRUE)
  e = me$all$eqtls
  e$UniqueID = paste0(e$snps, ".", e$gene)
  e = e[match(eqtl$UniqueID, e$UniqueID),]
  rownames(e) = e$UniqueID
  e[,c(3,4,6)]
}, mc.cores=1) #mc cores corresponds to number of cores available for multicore processing
gStats = do.call("cbind", gStatList) #combining all these different subgroup analyses into one
colnames(gStats) = gsub("\\.", "_", colnames(gStats)) #replacing periods with underscores
eqtl = cbind(eqtl, gStats) #add the subgroup eQTL analysis
colnames(eqtl)[-(1:3)] = paste0("Linear_" , colnames(eqtl)[-(1:3)])
#### full eqtl run ANOVA ####
 meQtl_all_anova = Matrix_eQTL_main(snps=snpSlice, 
	 gene = exprs, 
	 cvrt = covs, 
	 output_file_name.cis =  "cis.txt" ,
	 output_file_name =  "trans.txt" ,
	 pvOutputThreshold.cis = 0,
	 pvOutputThreshold=1,
	 snpspos = snpspos,
	 genepos = posExprs,
	 useModel = modelANOVA,
	 cisDist=1e6,
	 pvalue.hist = 100,
	 min.pv.by.genesnp = TRUE)
	 
eqtl_anova = meQtl_all_anova$all$eqtls
eqtl_anova$UniqueID = paste0(eqtl_anova$snps, ".", eqtl_anova$gene)
eqtl_anova = eqtl_anova[match(eqtl$UniqueID, eqtl_anova$UniqueID),]
rownames(eqtl_anova) = eqtl_anova$UniqueID
eqtl_anova = eqtl_anova[,3:4]
colnames(eqtl_anova) = paste0( "All_", colnames(eqtl_anova)   )
###########################################
########## anova subgroups ###################
###########################################

#### Splitting Into Subgroups ####

#### eQTL Analysis for Each Subgroup Defined Above
#Now we just apply a loop to perform eQTL analysis for each group defined above. 
gStatList = parallel::mclapply(gIndexes, function(ii) {
  yy = SlicedData$new( as.matrix(yExprs[,ii]) )
  ssnp = SlicedData$new( as.matrix(snpAll[,ii]) )
  covs = SlicedData$new(t(mod[ii,-1]))
  me = Matrix_eQTL_main(snps=ssnp, 
                        gene = yy, 
                        cvrt = covs,
                        output_file_name.cis =  "cis.txt" ,
						output_file_name =  "trans.txt" ,
                        pvOutputThreshold.cis = 0,  
                        pvOutputThreshold=1,
                        snpspos = snpspos, 
                        genepos = posExprs, 
                        useModel = modelANOVA,
						cisDist=1e6,
                        pvalue.hist = 100,
                        min.pv.by.genesnp = TRUE)
  e = me$all$eqtls
  e$UniqueID = paste0(e$snps, ".", e$gene)
  e = e[match(eqtl$UniqueID, e$UniqueID),]
  rownames(e) = e$UniqueID
  e[,c(3,4)]
}, mc.cores=1) #mc cores corresponds to number of cores available for multicore processing
gStats = do.call("cbind", gStatList) #combining all these different subgroup analyses into one
colnames(gStats) = gsub("\\.", "_", colnames(gStats)) #replacing periods with underscores

###
eqtl_anova = cbind(eqtl_anova,gStats)
colnames(eqtl_anova)[] = paste0("ANOVA_" , colnames(eqtl_anova)[])

### Combining ANOVA and Linear model results and saving
eqtl = cbind(eqtl, eqtl_anova) #add the subgroup eQTL analysis
save(eqtl, file = '/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Caudate/rdas/rawEqtls_Caudate_RiboZero_age13_merged.rda')