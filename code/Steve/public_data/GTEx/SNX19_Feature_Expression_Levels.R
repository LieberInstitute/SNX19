#SNX19 expression levels
library(GenomicRanges)
library(matrixStats)
library(MatrixEQTL)
library(jaffelab)
library(parallel)

foi <- read.csv('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_10Features_of_Interest_crossBuild_03_22_18.csv',stringsAsFactors=FALSE)
##########################
### load GTEx data
load('/users/ajaffe/Lieber/Projects/RNAseq/SNX19/GTEX/rawAndRpkmCounts_plusGenotype_SNX19_GTEX_n8206.rda')

#### Filtering out lowly expressed features ####
gIndex=which(rowMeans(geneRpkm) > 0.01) #selecting out genes with sufficient expression (genes with less than 0.01 read per kilobase million are not perserved)
geneRpkm = geneRpkm[gIndex,] #keeping only the genes above the rpkm threshold (threshold 0.01 here)
geneMap = geneMap[gIndex,] #keeping only the gene map for genes above the rpkm threshold
#geneMap is annotation information on the genes present (ensembl gene id, chromosome number, start site, end site, strand (+/-), length, HGNC symbol, and Entrez ID)
eIndex=which(rowMeans(exonRpkm) > 0.01)  #selecting out exons with sufficient expression (exons with less than 0.01 read per kilobase million are not perserved)
exonRpkm = exonRpkm[eIndex,] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
exonMap = exonMap[eIndex,] #keeping only the gene map for genes above the rpkm threshold
#geneMap is annotation information on the genes present (ensembl gene id, chromosome number, start site, end site, strand (+/-), length, HGNC symbol, and Entrez ID)
jIndex = which(rowMeans(jRpkm) > 0.2)#There are an excessive number of junctions. Apply a stricter filter and drop novel junctions to reduce the number.
jRpkm = jRpkm[jIndex,] #keeping only the junctions above the rpkm threshold (threshold = 0.01 rpkm here) and are NOT novel
jMap= jMap[jIndex] #keeping only the junction map for junctions above the rpkm threshold and are NOT novel


#### Filtering down to CAUC only ####
pIndexes <- pd$RACE==3
pd= pd[pIndexes,] #Keeping phenotypic data information only for the subjects above the age of 13
snp2 = as.matrix(snpAll[,pIndexes]) #Keeping SNP information only for subjects above the age of 13
geneRpkm2 = as.matrix(log2(geneRpkm[,pIndexes]+1)) #scaling gene rpkm data
exonRpkm2 = as.matrix(log2(exonRpkm[,pIndexes]+1)) #scaling exon rpkm data
jRpkm2 = as.matrix(log2(jRpkm[,pIndexes]+1)) #scaling junction rpkm data
exprs3features = rbind(geneRpkm2, exonRpkm2, jRpkm2) #combine the three types of data into a single matrix

exprs_foi<- exprs3features[rownames(exprs3features)%in%foi$Feature,]
pd2 = cbind(pd,t(exprs_foi))
pd2 = tidyr::gather(pd2, Feature, Normalized_Counts, (ncol(pd2)-nrow(foi)+1):ncol(pd2))
pd2 = cbind(pd2 , Tag=foi$Tag[match(pd2$Feature, foi$Feature)] )


feature_sum <- dplyr::group_by(pd2, SMTSD, Feature)
feature_sum$Normalized_Counts = feature_sum$Normalized_Counts
feature_sum_brain = feature_sum[grepl("Brain",feature_sum$SMTSD), ]
#plotting
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 axis.text=element_text(size=22),
				 axis.title=element_text(size=22),
				 legend.position="none") ) 
#a = ggplot(data=feature_sum_brain,aes(x=SMTSD,y=Normalized_Counts,fill=SMTSD)) + geom_boxplot() + facet_wrap(~Feature,scales='free')
#ggsave(a,filename='/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/plots/brain_8feature_distribution.pdf')
library(RColorBrewer)
n <- 1
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

##
#fix the labels for
feature_sum_brain$BrainRegion = gsub("Brain - ", "", feature_sum_brain$SMTSD)
##
ordered_feature_tag = c('junc8.10',
'junc8.8a',
'junc8c.9',
'ENSE00001757032.2',
'ENSE00002148491.1',
'junc1.2_long',
'junc2.3_short',
'junc4.6',
'junc8.9',
'SNX19 Gene')
ordered_feature = foi[match(foi$Tag, ordered_feature_tag),'Feature']
 pdf('/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/plots/brain_8feature_distribution.pdf',width=11,height=8.5)
	 for (k in ordered_feature ) {
	 dat2 = subset(feature_sum_brain, Feature==k)
  p <- ggplot(dat2, aes(x=BrainRegion, y=Normalized_Counts, col=BrainRegion, fill=BrainRegion)) + 
	   geom_boxplot(outlier.colour = NA, alpha = 0.3, aes(fill = BrainRegion)) + 
		geom_point(alpha = 1, aes(colour = BrainRegion), position = position_jitterdodge(jitter.width = NULL, jitter.height = 0,
  dodge.width = 0.75)
) + 
		labs(x = "Brain Region", 
			 y = "log2(Normalized Expression+1)",
			 title = paste0(unique(dat2$Tag),"\n",unique(dat2$Feature)) ) +
		theme(legend.position="none",axis.title.x=element_blank(), axis.text.x  = element_text(angle=45,size=14,hjust=1), axis.text.y =element_text(size=14) )
  print(p)
}
 dev.off()
 
feature_sum$SMTSD = as.factor(feature_sum$SMTSD)
x_axis_col <- ifelse(grepl("Brain",levels(feature_sum$SMTSD)), "blue", "black")

  pdf('/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/plots/49tissue_8feature_distribution.pdf',width=11,height=8.5)
	 for (k in ordered_feature ) {
	 dat2 = subset(feature_sum, Feature==k)
  p <- ggplot(dat2, aes(x=SMTSD, y=Normalized_Counts, col=SMTSD, fill=SMTSD)) + 
	   geom_boxplot(outlier.colour = NA, alpha = 0.3, aes(fill = SMTSD)) + 
		labs(x = "Brain Region", 
			 y = "log2(Normalized Expression+1)",
			 title = paste0(unique(dat2$Tag),"\n",unique(dat2$Feature)) ) + 
		theme(legend.position="none",axis.title.x=element_blank(), axis.text.x  = element_text(angle=60,size=10,hjust=1), axis.text.y =element_text(size=18) ) +
    theme(axis.text.x=element_text(colour = x_axis_col))
  print(p)
}
 dev.off()