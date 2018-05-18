#Load libaries
library('GenomicRanges')
library(AnnotationHub)
library(rtracklayer)

theme_set(theme_bw(base_size=14) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="bottom"	)) 
#load features of interest
foi = read.csv('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest.csv',stringsAsFactors=FALSE)
foi$Feature_hg19 = gsub("\\(-\\)", "(*)", foi$Feature)


### Cortecon hg19 analysis
load('/users/ajaffe/Lieber/Projects/RNAseq/LibdStem/rdas/junctionCounts.rda')
load('/users/ajaffe/Lieber/Projects/RNAseq/LibdStem/rdas/annotated_pheno_LIBD_stemCell.rda')
jCounts = juncCounts$countDF
jMap = juncCounts$anno
ll = which(rownames(jCounts) %in%foi$Feature_hg19)
snx19_jcounts = as.data.frame(jCounts[ll,])
snx19_jmap = as.data.frame(jMap[ll,])
snx19_jmap$Tag = foi[match(rownames(snx19_jmap),foi$Feature_hg19 ), 'Tag']
#juncCounts
dat = cbind(pd, t(snx19_jcounts) )
dat$Day = as.factor(dat$Day)
pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/eight_feature/cortecon.pdf')
for (i in rownames(snx19_jcounts)  ) {
dat2 = dat[,c('Day','Treatment',i)]
colnames(dat2)[3] <- "Junction"
a = ggplot(dat2, aes(x = Day, y = Junction)) + labs(x="Day",y="Junction Counts", title = i) + geom_boxplot(aes(fill=Treatment),outlier.colour = NA, alpha = .4 ) + 
		geom_point(aes(col=Treatment),alpha = 1, position=position_jitterdodge(jitter.width =0.2,dodge.width = 0.75) )
print(a)
}
dev.off()

#Build mapping table hg19 SNX19 to hg38 SNX19
load('/dcl01/lieber/RNAseq/23andMe/23andMe_analysis/eQTL_runs/Feature_Map_hg38_to_hg19.rda')
feature_lift$Hg19_FeatureName = gsub("\\..*","",feature_lift$Hg19_FeatureName)
feature_lift$Hg38_FeatureName = gsub("\\..*","",feature_lift$Hg38_FeatureName)
feature_lift[feature_lift$Type=="Gene", 'Shared_FeatureName'] = gsub("\\..*","",feature_lift[feature_lift$Type=="Gene", 'Shared_FeatureName'])



feature_lift[match(foi$Feature_hg19, feature_lift$Hg19_FeatureName),]

kk = match(foi$Feature, feature_lift$Hg19_FeatureName)
feature_lift[kk,]
foi[is.na(kk),]

#
load('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/junction_txdb/junction_annotation_hg19_gencode_v25lift37.rda')
hg19_jxn = theJunctions
names(hg19_jxn) = paste0(seqnames(hg19_jxn), ":", start(hg19_jxn), "-", end(hg19_jxn ), "(",  "*", ")" ) 
SNX19_hg19_gencode_v25lift37jxn = as.data.frame(hg19_jxn[hg19_jxn$symbol=="SNX19"])
#SNX19_hg19_gencode_v25lift37jxn$tx = sapply(SNX19_hg19_gencode_v25lift37jxn$tx, function(x) paste(x, sep=';') )

write.csv(SNX19_hg19_gencode_v25lift37jxn[,!colnames(SNX19_hg19_gencode_v25lift37jxn) %in% 'tx'],'/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_hg19_gencode_v25lift37jxn.csv')

foi[!foi$Feature_hg19 %in% names(hg19_jxn),]

### Lifting the junctions by hand


### Zhang: cell Type
1) Zhang: /dcl01/lieber/ajaffe/PublicData/Brain/Zhang_CellType (http://www.cell.com/neuron/abstract/S0896-6273(15)01019-3)


