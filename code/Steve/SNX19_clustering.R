## Load features of interest
load('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild.rda')

#### eQTL heatmaps

load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417.rda')
## Top 10 SNPs per region per transcript feature -duplicates =>  (400-dup)
library(tidyr)
library(dplyr)
tmp = cisMerge[,1:203]
tmp = tmp[,-13]
tmp = tmp[,-c(49:53)]
##### Tidying up the stats object
tidyMergedStats=  tmp%>%
  gather(key, value, 13:ncol(tmp))

###
tidyMergedStats$Model <- gsub("\\_[^\\_]*$","",tidyMergedStats[,'key'], )
tidyMergedStats$Statistic = gsub("(^.*_)","",tidyMergedStats[,'key'], )
###
tidyStats = tidyMergedStats[,-grep("key",colnames(tidyMergedStats))]  
tidyStats = spread(tidyStats,Statistic,value)

###
tidyStats$Region =gsub("_Linear.*","",tidyStats$Model)
tidyStats$Region =gsub("_ANOVA.*","",tidyStats$Region)

tidyStats$model =  mapply(function(p, r, x) sub(p, r, x, fixed = TRUE), p=tidyStats$Region, r=rep("",length(tidyStats$Region)), x=tidyStats$Model) 
tidyStats$model = substring(tidyStats$model, 2)
#save(tidyStats,file='/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/tidyStats_FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417.rda')
###
bestHits = tidyStats %>% filter(model=="Linear_All"& Feature %in% foi$Feature) %>% group_by(Region,Feature) %>%  top_n(n = -10, wt = pvalue) %>% as.data.frame()
bestHits$UniqueID = paste0(bestHits$SNP, ".", bestHits$Feature)
###
tStats = cisMerge[cisMerge$UniqueID%in% unique(bestHits$UniqueID),grep("Linear_All_statistic",colnames(cisMerge),value=T)]
rownames(tStats) =  cisMerge[cisMerge$UniqueID%in% unique(bestHits$UniqueID), 'UniqueID']

annotation <- data.frame( TranscriptFeature = jaffelab::ss(rownames(tStats),"\\.",2))
rownames(annotation) <- rownames(tStats)

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/clusteredEuclidean_eqtl_tvalues_257uniqueEqtls_fromBest10SNPs.pdf',height=20,width=60,onefile=F)
pheatmap::pheatmap(t(tStats),annotation_col = annotation, clustering_distance_cols = "euclidean", fontsize = 14)
dev.off()

######## Gwas sig best hits
bestHitsGwasSig = tidyStats %>% filter(model=="Linear_All"& Feature %in% foi$Feature & !is.na(pgc52_P) & pgc52_P<5e-8 ) %>% group_by(Region,Feature) %>%  top_n(n = -10, wt = pvalue) %>% as.data.frame()

bestHitsGwasSig$UniqueID = paste0(bestHitsGwasSig$SNP, ".", bestHitsGwasSig$Feature)
###
tStats = cisMerge[cisMerge$UniqueID%in% unique(bestHitsGwasSig$UniqueID),grep("Linear_All_statistic",colnames(cisMerge),value=T)]
rownames(tStats) =  cisMerge[cisMerge$UniqueID%in% unique(bestHitsGwasSig$UniqueID), 'UniqueID']

annotation <- data.frame( TranscriptFeature = jaffelab::ss(rownames(tStats),"\\.",2))
rownames(annotation) <- rownames(tStats)

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/clusteredEuclidean_eqtl_tvalues_257uniqueEqtls_fromBest10SNPs_gwasSigOnly.pdf',height=20,width=60,onefile=F)
pheatmap::pheatmap(t(tStats),annotation_col = annotation, clustering_distance_cols = "euclidean", fontsize = 14)
dev.off()
openxlsx::write.xlsx(list(anySNP=bestHits, gwasSig = bestHitsGwasSig ), file='/dcl01/lieber/ajaffe/Steve/SNX19/tables/Best_Eqtls_byRegion_byFeature_top10_tiesIncluded.xlsx' )

###### Another way of doing the heatmaps
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/tidyStats_FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC_100417.rda')
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/Annotated_5LIBD-CMC_set_Merged.rda')
library(dplyr)
library(tidyr)

###
bestHits = tidyStats %>% filter(model=="Linear_All"& Feature %in% foi$Feature) %>% 
						 group_by(Region,Feature) %>%  
						 top_n(n = -10, wt = pvalue) %>% 
						 as.data.frame()				 
bestHits$UniqueID = paste0(bestHits$SNP, ".", bestHits$Feature)

## Subsetting big object
tidyCisMerge = cisMerge[cisMerge$SNP %in% bestHits$SNP & cisMerge$Feature %in% foi$Feature, ]
tidyCisMerge = tidyCisMerge[,c('SNP',"Feature",grep("Linear_All_statistic",colnames(cisMerge),value=T))]
colnames(tidyCisMerge) = gsub("_Linear.*","",colnames(tidyCisMerge))
tidyCisMerge = gather(tidyCisMerge, Region,Tstat, -SNP,-Feature)
tidyCisMerge$Feature <- foi[match(tidyCisMerge$Feature, foi$Feature),'Tag']				
tidyCisMerge = unite(tidyCisMerge, Region_Feature, Region,Feature,sep="-")
tidyCisMerge = spread(tidyCisMerge, Region_Feature,Tstat)
###
tStats = tidyCisMerge
rownames(tStats) =  tStats$SNP
tStats = tStats[,-1]

annotation <- data.frame( Feature = jaffelab::ss(colnames(tStats),"-",2),
						  Dataset = ifelse(grepl("CMC",colnames(tStats)),"CMC","LIBD"),
						  Library = ifelse(grepl("PolyA",colnames(tStats)),"PolyA+","RiboZero"),
						  Region = jaffelab::ss(colnames(tStats),"-",1),
						  stringsAsFactors=F)
annotation$Region[grep("DLPFC",annotation$Region)] <- "DLPFC"				  
#annotation$Feature <- foi[match(annotation$Feature, foi$Feature),'Tag']		
annotation <- as.data.frame(unclass(annotation))					 				  
rownames(annotation) <- colnames(tStats)
annotation=annotation[,c('Feature','Region','Dataset','Library')]
library(RColorBrewer)
ann_colors = list(Feature = brewer.pal(8,"Dark2"),
				  Dataset = c('black','grey'),
				  Library = brewer.pal(5,"Set1")[4:5],
				  Region = brewer.pal(5,"Set1")[1:3] )
names(ann_colors[['Feature']])<- levels(annotation$Feature)
names(ann_colors[['Dataset']])<- levels(annotation$Dataset)
names(ann_colors[['Library']])<- levels(annotation$Library)
names(ann_colors[['Region']])<-levels(annotation$Region)
	
breaksList = seq(-1, 1, length.out =100)

sort_hclust <- function(...) as.hclust(dendsort::dendsort(as.dendrogram(...)))
mat_cluster_cols <- hclust(as.dist( 1- cor(tStats,use="pairwise.complete.obs") ) )
mat_cluster_cols <- sort_hclust(mat_cluster_cols)

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/clusteredCorrelation_eqtl_tvalues_179uniqueSNPs_fromBest10SNPs_regionFeatureTogether.pdf',height=20,width=30,onefile=T)
plot(mat_cluster_cols,cex=2.5)	

par(mar=c(30, 4, 4, 2) + 0.1	)
rafalib::myplclust(mat_cluster_cols, lab=mat_cluster_cols$labels, lab.col=ann_colors$Feature[as.character(annotation$Feature)],hang=.1,cex=2.5 )
par(mar=c(30, 4, 4, 2) + 0.1	)
rafalib::myplclust(mat_cluster_cols, lab=mat_cluster_cols$labels, lab.col=ann_colors$Dataset[as.character(annotation$Dataset)],hang=.1,cex=2.5)
par(mar=c(30, 4, 4, 2) + 0.1	)
rafalib::myplclust(mat_cluster_cols, lab=mat_cluster_cols$labels, lab.col=ann_colors$Library[as.character(annotation$Library)],hang=.1,cex=2.5 )
par(mar=c(30, 4, 4, 2) + 0.1	)
rafalib::myplclust(mat_cluster_cols, lab=mat_cluster_cols$labels, lab.col=ann_colors$Region[as.character(annotation$Region)],hang=.1,cex=2.5)

pheatmap::pheatmap(cor(tStats,use="pairwise.complete.obs"),
				   breaks=breaksList,
				   clustering_distance_cols = "euclidean",
				   clustering_distance_rows = "euclidean",
				   cluster_cols= mat_cluster_cols,
				   cluster_rows= mat_cluster_cols,
				   treeheight_row=200,
				   treeheight_col=0,
				   annotation_row = annotation, 				   
				   annotation_colors = ann_colors, 
				   show_colnames = FALSE, 
				   show_rownames = TRUE,
				   fontsize = 14)
dev.off()	 
	 
###	 
pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/clusteredEuclidean_eqtl_tvalues_179uniqueSNPs_fromBest10SNPs_regionFeatureTogether.pdf',height=20,width=60,onefile=F)
pheatmap::pheatmap(t(tStats),annotation_row = annotation, clustering_distance_cols = "euclidean", fontsize = 14)
dev.off()

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/clusteredCorrelation_eqtl_tvalues_179uniqueSNPs_fromBest10SNPs_SNPcor.pdf',height=40,width=40,onefile=F)
pheatmap::pheatmap(cor(t(tStats),use="pairwise.complete.obs"), clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", fontsize = 14)
dev.off()


######### GWAS Sig

bestHits_correct = tidyStats %>% filter(model=="Linear_All"& Feature %in% foi$Feature & !is.na(pgc52_P) & pgc52_P<5e-8 ) %>% 
						 group_by(Region,Feature) %>%  
						 top_n(n = -10, wt = pvalue) %>% 
						 as.data.frame()				 

bestHits_incorrect = tidyStats %>% filter(model=="Linear_All"& Feature %in% foi$Feature & !is.na(pgc52_P) & pgc52_P<5e-8 ) %>% 
						 group_by(Region,Feature) %>%  
						 top_n(n = 10, wt = pvalue) %>% 
						 as.data.frame()	
gwas_snp = tidyStats %>% filter(model=="Linear_All"& Feature %in% foi$Feature & !is.na(pgc52_P) & pgc52_P<5e-8 )
library(eulerr)
fit <- euler(list( CORRECT_SNPS=unique(bestHits_correct$SNP), INCORRECT_SNPS = unique(bestHits_incorrect$SNP), GWAS_SNPS_TESTED=unique(gwas_snp$SNP) )	)

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/venn_diagram_of_correct_vs_incorrect_heatmap_SNPs.pdf')
plot(fit,quantities=T, fill_opacity = 0.3, fill=c('red','blue','green') )
dev.off()		 
## Subsetting big object
tidyCisMerge = cisMerge[cisMerge$SNP %in% bestHits$SNP & cisMerge$Feature %in% foi$Feature, ]
tidyCisMerge = tidyCisMerge[,c('SNP',"Feature",grep("Linear_All_statistic",colnames(cisMerge),value=T))]
colnames(tidyCisMerge) = gsub("_Linear.*","",colnames(tidyCisMerge))
tidyCisMerge = gather(tidyCisMerge, Region,Tstat, -SNP,-Feature)
tidyCisMerge$Feature <- foi[match(tidyCisMerge$Feature, foi$Feature),'Tag']				
tidyCisMerge = unite(tidyCisMerge, Region_Feature, Region,Feature,sep="-")
tidyCisMerge = spread(tidyCisMerge, Region_Feature,Tstat)
###
tStats = tidyCisMerge
rownames(tStats) =  tStats$SNP
tStats = tStats[,-1]

annotation <- data.frame( Feature = jaffelab::ss(colnames(tStats),"-",2),
						  Dataset = ifelse(grepl("CMC",colnames(tStats)),"CMC","LIBD"),
						  Library = ifelse(grepl("PolyA",colnames(tStats)),"PolyA+","RiboZero"),
						  Region = jaffelab::ss(colnames(tStats),"-",1),
						  stringsAsFactors=F)
annotation$Region[grep("DLPFC",annotation$Region)] <- "DLPFC"				  
#annotation$Feature <- foi[match(annotation$Feature, foi$Feature),'Tag']		
annotation <- as.data.frame(unclass(annotation))					 				  
rownames(annotation) <- colnames(tStats)
annotation=annotation[,c('Feature','Region','Dataset','Library')]

library(RColorBrewer)
ann_colors = list(Feature = brewer.pal(8,"Dark2"),
				  Dataset = c('black','grey'),
				  Library = brewer.pal(5,"Set1")[4:5],
				  Region = brewer.pal(5,"Set1")[1:3] )
names(ann_colors[['Feature']])<- levels(annotation$Feature)
names(ann_colors[['Dataset']])<- levels(annotation$Dataset)
names(ann_colors[['Library']])<- levels(annotation$Library)
names(ann_colors[['Region']])<-levels(annotation$Region)
	
breaksList = seq(-1, 1, length.out =100)

sort_hclust <- function(...) as.hclust(dendsort::dendsort(as.dendrogram(...)))
mat_cluster_cols <- hclust(as.dist( 1- cor(tStats,use="pairwise.complete.obs") ) )
mat_cluster_cols <- sort_hclust(mat_cluster_cols)


pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/clusteredCorrelation_eqtl_tvalues_102uniqueSNPs_fromBest10SNPs_regionFeatureTogether_gwasSigSNP.pdf',height=20,width=30,onefile=T)
plot(mat_cluster_cols,cex=2.5)

par(mar=c(35, 4, 4, 2) + 0.1	)
rafalib::myplclust(mat_cluster_cols, lab=mat_cluster_cols$labels, lab.col=ann_colors$Feature[as.character(annotation$Feature)],hang=.1,cex=2.3 )
par(mar=c(35, 4, 4, 2) + 0.1	)
rafalib::myplclust(mat_cluster_cols, lab=mat_cluster_cols$labels, lab.col=ann_colors$Dataset[as.character(annotation$Dataset)],hang=.1,cex=2.3)
par(mar=c(35, 4, 4, 2) + 0.1	)
rafalib::myplclust(mat_cluster_cols, lab=mat_cluster_cols$labels, lab.col=ann_colors$Library[as.character(annotation$Library)],hang=.1,cex=2.3 )
par(mar=c(35, 4, 4, 2) + 0.1	)
rafalib::myplclust(mat_cluster_cols, lab=mat_cluster_cols$labels, lab.col=ann_colors$Region[as.character(annotation$Region)],hang=.1,cex=2.3)
			   
pheatmap::pheatmap(cor(tStats,use="pairwise.complete.obs"),
				   breaks=breaksList,
				   clustering_distance_cols = "euclidean",
				   clustering_distance_rows = "euclidean",
				   cluster_cols= mat_cluster_cols,
				   cluster_rows= mat_cluster_cols,
				   treeheight_row=200,
				   treeheight_col=0,
				   annotation_row = annotation, 				   
				   annotation_colors = ann_colors, 
				   show_colnames = FALSE, 
				   show_rownames = TRUE,
				   fontsize = 14)
dev.off()

######
pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/clusteredEuclidean_eqtl_tvalues_102uniqueSNPs_fromBest10SNPs_regionFeatureTogether_gwasSigSNP.pdf',height=20,width=60,onefile=F)
pheatmap::pheatmap(t(tStats),annotation_row = annotation, clustering_distance_cols = "euclidean", fontsize = 14, annotation_colors = ann_colors,)
dev.off()

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/clusteredCorrelation_eqtl_tvalues_102uniqueSNPs_fromBest10SNPs_SNPcor_gwasSigSNP.pdf',height=40,width=40,onefile=F)
pheatmap::pheatmap(cor(t(tStats),use="pairwise.complete.obs"), clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", fontsize = 14, breaks=breaksList)
dev.off()

############ Expression heatmaps

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

#######
shared_people = Reduce(intersect,list(pd_DlpfcPolyA$BrNum, pd_DlpfcRiboZero$BrNum, pd_Caudate$BrNum, pd_Hippo$BrNum) )

##
cut_DlpfcPolya = yExprs_DlpfcPolya[foi$Feature, pd_DlpfcPolyA[match(shared_people, pd_DlpfcPolyA$BrNum), 'RNum'] ]
colnames(cut_DlpfcPolya) <- paste0("DlpfcPolya_",shared_people )
##
cut_DlpfcRiboZero = yExprs_DlpfcRiboZero[foi$Feature, pd_DlpfcRiboZero[match(shared_people, pd_DlpfcRiboZero$BrNum), 'RNum'] ]
colnames(cut_DlpfcRiboZero) <- paste0("DlpfcRiboZero_",shared_people )
##
cut_Caudate = yExprs_Caudate[foi$Feature, pd_Caudate[match(shared_people, pd_Caudate$BrNum), 'RNum'] ]
colnames(cut_Caudate) <- paste0("Caudate_",shared_people )
##
cut_Hippo = yExprs_Hippo[foi$Feature, pd_Hippo[match(shared_people, pd_Hippo$BrNum), 'RNum'] ]
colnames(cut_Hippo) <- paste0("Hippo_",shared_people )
##
allRegions = cbind(cut_DlpfcPolya, cut_DlpfcRiboZero, cut_Caudate,cut_Hippo )


annotation <- data.frame( Region = jaffelab::ss(colnames(allRegions),"_",1),
						  Dx=pd[match(shared_people, pd$BrNum ),"Dx"], 
						  Sex = pd[match(shared_people, pd$BrNum ),"Sex"],
						  Race = pd[match(shared_people, pd$BrNum ),"Race"] )
rownames(annotation) <- colnames(allRegions)

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/clusteredCorrelation_log2Expression_n217_adults_sharedBrains.pdf',height=20,width=60,onefile=F)
pheatmap::pheatmap(allRegions,annotation_col = annotation, labels_col = jaffelab::ss(colnames(allRegions),"_",2), clustering_distance_cols = "correlation", fontsize_row = 20, fontsize_col = 5,scale='row')
dev.off()

###
foi_DlpfcPolya = rowMeans(yExprs_DlpfcPolya[foi$Feature, ] )
foi_DlpfcRiboZero = rowMeans(yExprs_DlpfcRiboZero[foi$Feature, ] )
foi_Caudate = rowMeans(yExprs_Hippo[foi$Feature, ] )
foi_Hippo = rowMeans(yExprs_Caudate[foi$Feature, ] )
foi_CMC = rowMeans(yExprs_CMC[foi$Feature, ] )

foi_meanExprs = cbind(foi_DlpfcPolya, foi_DlpfcRiboZero, foi_Caudate,foi_Hippo,foi_CMC  )

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/clusteredCorrelation_mean_log2Expression_adults_samplesUsed_featureNormalized.pdf',height=10,width=10,onefile=F)
pheatmap::pheatmap(foi_meanExprs, labels_col = jaffelab::ss(colnames(foi_meanExprs),"_",2), labels_row = foi[match(foi$Feature,rownames(foi_meanExprs)),'Tag'], clustering_distance_cols = "correlation", clustering_distance_rows='correlation', fontsize_row = 15, fontsize_col = 15,scale='row')
dev.off()

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/clusteredCorrelation_mean_log2Expression_adults_samplesUsed_DatasetNormalized.pdf',height=10,width=10,onefile=F)
pheatmap::pheatmap(foi_meanExprs, labels_col = jaffelab::ss(colnames(foi_meanExprs),"_",2), labels_row = foi[match(foi$Feature,rownames(foi_meanExprs)),'Tag'], clustering_distance_cols = "correlation", clustering_distance_rows='correlation', fontsize_row = 15, fontsize_col = 15,scale='column')
dev.off()

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/clusteredCorrelation_mean_log2Expression_adults_samplesUsed_NotNormalized.pdf',height=10,width=10,onefile=F)
pheatmap::pheatmap(foi_meanExprs, labels_col = jaffelab::ss(colnames(foi_meanExprs),"_",2), labels_row = foi[match(foi$Feature,rownames(foi_meanExprs)),'Tag'], clustering_distance_cols = "correlation", clustering_distance_rows='correlation', fontsize_row = 15, fontsize_col = 15,scale='none')
dev.off()

############## sanity checks

##
x = tStats[,c(grep("SNX19.gene",colnames(tStats)),grep("junc8.9",colnames(tStats))) ]
x = x[,gtools::mixedsort(colnames(x))]
cor(x, use="pairwise.complete.obs")


##
library(ggplot2)
library(GGally)
x2 = spread(tidyCisMerge,Feature,Tstat)
pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/GGally_pairwise_scatter_8features.pdf',height=30,width=30)
ggplot(data=x2, aes(x=SNX19.gene, y=junc8.9) ) + geom_point() + facet_wrap(~Region)
GGally::ggscatmat(x2, columns = 3:10, color="Region", alpha=0.8)
dev.off()