library('GenomicRanges')
load('/dcl01/lieber/ajaffe/PublicData/Brain/Zhang_CellType/hg38/rpkmCounts_Zhang_CellType_Human_hg38_n41.rda')
load('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild.rda')
pd=data.table::fread('/dcl01/lieber/ajaffe/PublicData/Brain/Zhang_CellType/SraRunTable.txt')

#Convert hg38 eIDs to ENSEMBL exon IDs
load('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/exon_names_hg19_hg38.rda')
mm = match(rownames(exonMap),hg38_exons$hg38_eID )
rownames(exonRpkm) <- hg38_exons[mm ,'gencode_exonID' ]
rownames(exonMap) <- hg38_exons[mm ,'gencode_exonID' ]
## Fix gene rownames
rownames(geneRpkm) = gsub("\\..*","",rownames(geneRpkm) )
geneMap$gencodeID =  gsub("\\..*","",geneMap$gencodeID )
##
FOI = foi$MergeName
foi[foi$exprsType=="Junction",'Tag'][!foi[foi$exprsType=="Junction",'MergeName']%in%rownames(jRpkm)] #features not mapped

### get expression subset
geneRpkm =  geneRpkm[rownames(geneRpkm)%in%FOI, ,drop=FALSE]
exonRpkm =  exonRpkm[rownames(exonRpkm)%in%FOI, ]
jRpkm = jRpkm[rownames(jRpkm)%in%FOI, ]

### get map subset
geneMap =  geneMap[geneMap$gencodeID%in% FOI, ,drop=FALSE]
exonMap =  exonMap[rownames(exonMap)%in% FOI, ]
jMap =  as.data.frame(jMap[names(jMap)%in% FOI])

## merge jMap
jMap$newName = paste0(jMap$seqnames, ":", jMap$start, "-", jMap$end, "(*)")
jMap = jMap[!duplicated(jMap$newName),]

## common expression map to annotate
geneMap$EnsemblGeneID = rownames(geneMap)
geneMap$Class = "InEns"
geneMap$Type = "Gene"

colnames(exonMap)[7] = "EnsemblGeneID"
exonMap$Class = "InEns"
exonMap$Type = "Exon"


jMap$Type = "Junction"
colnames(jMap)[10] <- "EnsemblGeneID"
name = c("EnsemblGeneID", "Symbol", "Type", "Class")
exprsMap = rbind(as.data.frame(geneMap)[, name],
	as.data.frame(exonMap)[, name],
	as.data.frame(jMap)[, name])
###
all_3_features = rbind(geneRpkm, exonRpkm, jRpkm)
all_3_feature_map = exprsMap

#Combining and gathering these objects
all_3_features_gathered <- as.data.frame(all_3_features)
all_3_features_gathered$Rowname <- rownames(all_3_features_gathered)
all_3_features_gathered <- tidyr::gather(all_3_features_gathered, Sample, Rpkm, -Rowname) 
all_3_features_gathered <- cbind(all_3_features_gathered, pd[match(all_3_features_gathered$Sample,pd$Run_s),] )
all_3_features_gathered$Type <- ifelse(grepl("ENSG",all_3_features_gathered$Rowname), "Gene", ifelse(grepl(":",all_3_features_gathered$Rowname),"Jxn", "Exon") )
all_3_features_gathered <- cbind(all_3_features_gathered,all_3_feature_map[match(all_3_features_gathered$Rowname,rownames(all_3_feature_map) ),] )
####
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none")) 
#all_3_features_gathered$Stage = factor(all_3_features_gathered$Stage,levels=c("Pluripotency","Neural Differentation","Cortical Specification","Deep Layer Formation","Upper Layer Formation" ) )
all_3_features_gathered = subset(all_3_features_gathered,source_name_s!="")

 pdf('/dcl01/lieber/ajaffe/Steve/SNX19/zhang/plots/SNX19_8feature_zhang_cellType.pdf')
	 for (k in unique(all_3_features_gathered$Rowname) ) {
	 dat2 = subset(all_3_features_gathered, Rowname==k)
	  p <- ggplot(dat2, aes(x=cell_type_s, y=log2(Rpkm+1)) ) +
			geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_jitter(alpha = 1, position = position_jitter(width = 0.2, height = 0), aes(col=source_name_s ) ) +
			labs(x = "Cell Type", 
				 y = "log2(RPKM+1)",
				 title = foi[match(k,foi$MergeName ),'Tag']) + theme(legend.position="right") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") +theme(axis.text.x  = element_text(angle=45,size=10,hjust=1))
	  print(p)
	  }
 dev.off()
 
 