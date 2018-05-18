## Absolute isoform 
library(jaffelab)
library(MatrixEQTL)
library(SummarizedExperiment)
library(matrixStats)

## Load stats
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417.rda')

## Load features of intereset
load('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild.rda')

## DLPFC
setwd('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/DLPFC_polyA')
## load counts
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_polyA/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_n495_updateMap_03082017.rda")

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

## get expression PCs ####
exprs3features = log2(rbind(geneRpkm, exonRpkm, jRpkm)+1) #combine the three types of data into a single matrix
pcaexprs3features = prcomp(t(exprs3features+1)) #principle component decomposition on the transpose of that matrix + 1
all3PCs = pcaexprs3features$x[,1:15] #take the first 15 PCs from that matrix

## set up model matrix
pd$Dx = factor(pd$Dx, levels=c("Control","Schizo") )
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
yExprs = t(as.data.frame( rowSums( scale( t(exprs3features[foi$Feature[!foi$Tag%in%c('junc8.9','junc4.6','SNX19.gene')],] ) ) )) )
rownames(yExprs) <- "totalIsoformExpression"
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
	 output_file_name.cis =  "cis.txt",
	 output_file_name =  "trans.txt",
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

#
a = rowMeans(exprs3features[foi$Feature,])
names(a) <- foi$Tag
a


##### Boxplots
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))

dat = cbind(pd, t(exprs3features[foi$Feature,]),t(yExprs), t(snpAll) )
colnames(dat)<-plyr::mapvalues(colnames(dat),foi$Feature, foi$Tag)
dat = tidyr::gather_(dat,'transcriptFeature','expression', c(foi$Tag, 'totalIsoformExpression' ) )

plot_dat = dat[!dat$transcriptFeature%in%c("junc4.6","junc8.9","SNX19.gene"),]
colnames(plot_dat) =  plyr::mapvalues(colnames(dat),snpMapAll$SNP, snpMapAll$name)
plot_dat[,colnames(plot_dat)%in%snpMapAll$name] <- lapply(plot_dat[,colnames(plot_dat)%in%snpMapAll$name], factor) 


pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/eqtl_boxplots_isoform_panels.pdf',height=12,width=15)
for (i in eqtl[!duplicated(eqtl$All_beta),'SNP'][1:10]) {
#Change column name for ggplot2 to work
k=snpMapAll$name[match(i, snpMapAll$SNP)]
b = ggplot(data=plot_dat[!is.na(plot_dat[,k]),], aes_string(x =k, y = 'expression')) + 
		facet_wrap(~transcriptFeature,ncol=3,nrow=2,scales='free') +
        geom_boxplot(outlier.colour = NA, alpha = 0.1, col='black')  + 
		geom_jitter(height=0) + 
		labs(y="Normalized Expression", title = i)	
print(b)			}
dev.off()


#### Pca of 5 expression features
pcs = prcomp(t(exprs3features[foi$Feature[!foi$Tag%in%c('junc8.9','junc4.6','SNX19.gene')],] ))
getPcaVars(pcs)
pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/five_feature_pca_plots.pdf')
plot(pcs$x[,1],pcs$x[,2] )
plot(pcs$x[,2],pcs$x[,3] )
plot(pcs$x[,3],pcs$x[,4] )
plot(pcs$x[,4],pcs$x[,5] )

plot(pcs$x[,1],pcs$x[,2] )
plot(pcs$x[,2],pcs$x[,3] )
plot(pcs$x[,3],pcs$x[,4] )
plot(pcs$x[,4],pcs$x[,5] )
dev.off()

###### matrixEqtl of PC1

#### full eqtl run ####
 meQtl_all = Matrix_eQTL_main(snps=snpSlice, 
	 gene = SlicedData$new(as.matrix( t(pcs$x[,1]) )), 
	 cvrt = covs, 
	 output_file_name.cis =  "cis.txt",
	 output_file_name =  "trans.txt",
	 pvOutputThreshold.cis = 0,  
	 pvOutputThreshold=1,
	 snpspos = snpspos, 
	 genepos = posExprs, 
	 useModel = modelLINEAR,	
	 cisDist=1e6,
	 pvalue.hist = 100,
	 min.pv.by.genesnp = TRUE)
	 
eqtl2 = meQtl_all$all$eqtls #extracting cis eqtl information from the eQTL analysis
colnames(eqtl2)[1:2] = c("SNP","Feature") #adding a column called feature
eqtl2$Feature = as.character(eqtl2$Feature) #the feature column contains expression information
eqtl2$SNP = as.character(eqtl2$SNP) #converting the eqtl snps variable to a character variable
colnames(eqtl2)[colnames(eqtl2) %in% c('statistic','pvalue','FDR','beta')] = paste0("All_", colnames(eqtl2)[colnames(eqtl2) %in% c('statistic','pvalue','FDR','beta')]) #Adding "ALL_ to "statistics, p-value, FDR, and beta" to help distinguish from future analyses 
eqtl2$UniqueID = paste0(eqtl2$SNP, ".", eqtl2$Feature) #making a variable called Unique ID... Just a combination of SNP and feature
eqtl2 = eqtl2[,c('UniqueID', 'SNP', 'Feature', 'All_statistic', 'All_pvalue', 'All_FDR', 'All_beta')] #just reordering the columns of the dataframe	

cor.test(eqtl$All_beta, eqtl2$All_beta[match(eqtl$SNP, eqtl2$SNP )])
cor.test(eqtl$All_pvalue, eqtl2$All_pvalue[match(eqtl$SNP, eqtl2$SNP )])
cor.test(eqtl$All_statistic, eqtl2$All_statistic[match(eqtl$SNP, eqtl2$SNP )])

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/total_isoformQTL_v_PC1_5Feature_QTL.pdf')
plot(eqtl$All_beta, eqtl2$All_beta[match(eqtl$SNP, eqtl2$SNP )], xlab="Total Isoform QTL Beta", ylab="5 Feature PC1 QTL Beta")
plot(eqtl$All_statistic, eqtl2$All_statistic[match(eqtl$SNP, eqtl2$SNP )], xlab="Total Isoform QTL Statistic", ylab="5 Feature PC1 QTL Statistic")
plot(-log10(eqtl$All_pvalue), -log10( eqtl2$All_pvalue[match(eqtl$SNP, eqtl2$SNP )]), xlab="Total Isoform QTL -log10(P)", ylab="5 Feature PC1 QTL -log10(P)")
dev.off()

#######
best_SNP_by_feature = cisMerge[  cisMerge$Feature %in% foi$Feature[!foi$Tag%in%c('junc8.9','junc4.6','SNX19.gene')], ]
library(dplyr)
five_feature_bestSNP = group_by(best_SNP_by_feature, Feature ) %>% summarise( best_SNP = SNP[which.min(DLPFC_PolyA_Linear_All_pvalue)], p=min(DLPFC_PolyA_Linear_All_pvalue),statistic=DLPFC_PolyA_Linear_All_statistic[which.min(DLPFC_PolyA_Linear_All_pvalue)]   ) %>% as.data.frame()

five_feature_bestSNP = rbind(five_feature_bestSNP, 
					   setNames(eqtl[1,c('Feature','SNP','All_pvalue','All_statistic')], names(five_feature_bestSNP)),
					   setNames(eqtl2[1,c('Feature','SNP','All_pvalue','All_statistic')] , names(five_feature_bestSNP)) )
					   
five_feature_bestSNP$Feature[five_feature_bestSNP$Feature=="row1"] <- "Five Feature PC1"

stats_summary = best_SNP_by_feature[best_SNP_by_feature$SNP%in%five_feature_bestSNP$best_SNP,c('Feature','SNP','DLPFC_PolyA_Linear_All_pvalue', "DLPFC_PolyA_Linear_All_statistic")]
stats_summary = rbind(stats_summary, 
					   setNames(eqtl[eqtl$SNP %in% stats_summary$SNP,c('Feature','SNP','All_pvalue','All_statistic')], names(stats_summary)),
					   setNames(eqtl2[eqtl2$SNP %in% stats_summary$SNP,c('Feature','SNP','All_pvalue','All_statistic')] , names(stats_summary)) )
					   
library(tidyr)
pval_Stats = tidyr::spread(stats_summary[,c('Feature','SNP','DLPFC_PolyA_Linear_All_pvalue')], SNP, DLPFC_PolyA_Linear_All_pvalue)					   
pval_Stats$Feature[pval_Stats$Feature =="row1"] <- "Five Feature PC1"			

t_Stats = tidyr::spread(stats_summary[,c('Feature','SNP','DLPFC_PolyA_Linear_All_statistic')], SNP, DLPFC_PolyA_Linear_All_statistic)				
t_Stats$Feature[t_Stats$Feature =="row1"] <- "Five Feature PC1"		

## Get R^2
cor(t(snpAll[unique(five_feature_bestSNP$best_SNP),]), use = c("pairwise.complete.obs")) ^2

		   