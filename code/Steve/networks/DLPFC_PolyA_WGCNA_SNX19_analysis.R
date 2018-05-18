library('corrplot')

## Extract SNX19 features

# load dlpfc polya expression
load("/dcs01/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")

# load snx19 feature list
load('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild.rda', verbose=TRUE)
foi = foi[order(foi$exprsType),]
foi$Length = foi$exprsEnd - foi$exprsStart + 1

##
foi$Feature[!foi$Feature%in%rownames(exonRpkm) & foi$exprsType=="Exon"]				
eIndex309 = which(exonMap$Symbol=="SNX19" & exonMap$Length==309) ## exons match for length and tx
eIndex111 = which(exonMap$Symbol=="SNX19" & exonMap$Length==111)
foi = foi[order(foi$exprsType),]

##
snxRpkm = rbind(exonRpkm[c(eIndex309,eIndex111),,drop=F],
				geneRpkm[rownames(geneRpkm) %in% foi$Feature,,drop=F],
				jRpkm[rownames(jRpkm) %in% foi$Feature,,drop=F])
rownames(snxRpkm)[1:3] = foi$Feature[c(1:3)]
foi = foi[match(rownames(snxRpkm), foi$Feature),]

## Load in WGCNA modules gene-level
###### GENES ####

load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/wgcna/LIBD_modules_twoRuns.rda", verbose=T)
datMEAdj = netAdj$MEs
datMEAdj=datMEAdj[,gtools::mixedsort(colnames(datMEAdj))] # reorder columns 
datMEAdj=datMEAdj[,-1] ## remove unassigned

datMEQsva = netQsva$MEs
datMEQsva=datMEQsva[,gtools::mixedsort(colnames(datMEQsva))] # reorder columns 
datMEQsva_gene=datMEQsva[,-1] ## remove unassigned

###### EXONS ####
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/wgcna/LIBD_modules_exonLevel_qSVA.rda", verbose=T)
datMEQsva = netQsva$MEs
datMEQsva=datMEQsva[,gtools::mixedsort(colnames(datMEQsva))] # reorder columns 
datMEQsva_exon=datMEQsva[,-1] ## remove unassigned

###### JUNCTIONS ####
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/wgcna/LIBD_modules_junctionLevel_qSVA.rda", verbose=T)
datMEQsva = netQsva$MEs
datMEQsva=datMEQsva[,gtools::mixedsort(colnames(datMEQsva))] # reorder columns 
datMEQsva_jxn=datMEQsva[,-1] ## remove unassigned

###
load("/users/ajaffe/Lieber/Projects/RNAseq/SzControl_DE_paper/rdas/geneLevel_LIBD_qSVsAndMod.rda",verbose=T)
snxRpkm = snxRpkm[,rownames(mod)]
rownames(snxRpkm) <- foi[match(rownames(snxRpkm), foi$Feature),'Tag']

###### plots for adjusted wgcna gene-level
res = Hmisc::rcorr(as.matrix(cbind(datMEAdj,t(snxRpkm) )))
cormat_gene_Adj=res$r[colnames(datMEAdj),rownames(snxRpkm)]
pmat_gene_Adj = res$P[colnames(datMEAdj),rownames(snxRpkm)]

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/corrplot_SNX19_8Features_vs_DLPFC_PolyA_WGCNA_Adj_21Modules.pdf',height=11,width=8.5)
corrplot(cormat_gene_Adj, p.mat=pmat_gene_Adj, sig.level = .05)
pmat_gene_Adj[cormat_gene_Adj<0] =1
cormat_gene_Adj[cormat_gene_Adj<0]=0
corrplot(corr=cormat_gene_Adj, p.mat=pmat_gene_Adj, sig.level = .05)
dev.off()

####### plots for qsva wgcna gene-level
res = Hmisc::rcorr(as.matrix(cbind(datMEQsva_gene,t(snxRpkm) )))
cormat_gene_qsva=res$r[colnames(datMEQsva_gene),rownames(snxRpkm)]
pmat_gene_qsva = res$P[colnames(datMEQsva_gene),rownames(snxRpkm)]

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/corrplot_SNX19_8Features_vs_DLPFC_PolyA_WGCNA_Qsva_32Modules.pdf',height=11,width=8.5)
corrplot(cormat_gene_qsva, p.mat=pmat_gene_qsva, sig.level = .05)
pmat_gene_qsva[cormat_gene_qsva<0] =1
cormat_gene_qsva[cormat_gene_qsva<0]=0
corrplot(corr=cormat_gene_qsva, p.mat=pmat_gene_qsva, sig.level = .05)
dev.off()

####### plots for qsva wgcna exon-level
res = Hmisc::rcorr(as.matrix(cbind(datMEQsva_exon,t(snxRpkm) )))
cormat_exon_qsva=res$r[colnames(datMEQsva_exon),rownames(snxRpkm)]
pmat_exon_qsva = res$P[colnames(datMEQsva_exon),rownames(snxRpkm)]

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/corrplot_SNX19_8Features_vs_DLPFC_PolyA_WGCNA_Qsva_110Modules_ExonLevel.pdf',height=50,width=50)
corrplot(cormat_exon_qsva, p.mat=pmat_exon_qsva, sig.level = .05)
pmat_exon_qsva[cormat_exon_qsva<0] =1
cormat_exon_qsva[cormat_exon_qsva<0]=0
corrplot(corr=cormat_exon_qsva, p.mat=pmat_exon_qsva, sig.level = .05)
dev.off()

res = Hmisc::rcorr(as.matrix(cbind(datMEQsva_exon,t(snxRpkm) )))
cormat_exon_qsva=res$r[colnames(datMEQsva_exon),rownames(snxRpkm)]
pmat_exon_qsva = res$P[colnames(datMEQsva_exon),rownames(snxRpkm)]

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/corrplot_SNX19_8Features_vs_DLPFC_PolyA_WGCNA_Qsva_110Modules_ExonLevel_sigOnly.pdf',height=10,width=20)
corrplot(cormat_exon_qsva[rowSums(pmat_exon_qsva<0.05)>0,], p.mat=pmat_exon_qsva[rowSums(pmat_exon_qsva<0.05)>0,], sig.level = .05)
dev.off()

####### plots for qsva wgcna junction-level
res = Hmisc::rcorr(as.matrix(cbind(datMEQsva_jxn,t(snxRpkm) )))
cormat_jxn_qsva=res$r[colnames(datMEQsva_jxn),rownames(snxRpkm)]
pmat_jxn_qsva = res$P[colnames(datMEQsva_jxn),rownames(snxRpkm)]

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/corrplot_SNX19_8Features_vs_DLPFC_PolyA_WGCNA_Qsva_228Modules_JunctionLevel.pdf',height=50,width=50)
corrplot(cormat_jxn_qsva, p.mat=pmat_jxn_qsva, sig.level = .05)
pmat_jxn_qsva[cormat_jxn_qsva<0] =1
cormat_jxn_qsva[cormat_jxn_qsva<0]=0
corrplot(corr=cormat_jxn_qsva, p.mat=pmat_jxn_qsva, sig.level = .05)
dev.off()

res = Hmisc::rcorr(as.matrix(cbind(datMEQsva_jxn,t(snxRpkm) )))
cormat_jxn_qsva=res$r[colnames(datMEQsva_jxn),rownames(snxRpkm)]
pmat_jxn_qsva = res$P[colnames(datMEQsva_jxn),rownames(snxRpkm)]

pdf('/dcl01/lieber/ajaffe/Steve/SNX19/plots/corrplot_SNX19_8Features_vs_DLPFC_PolyA_WGCNA_Qsva_228Modules_JunctionLevel_sigOnly.pdf',height=30,width=50)
corrplot(cormat_jxn_qsva[rowSums(pmat_jxn_qsva<0.05)>0,], p.mat=pmat_jxn_qsva[rowSums(pmat_jxn_qsva<0.05)>0,], sig.level = .05)
dev.off()


### Save out correlation matrix for significant rows
openxlsx::write.xlsx(list(geneCor_Pearson_R_AdjNetwork = tibble::rownames_to_column( as.data.frame( cormat_gene_Adj[rowSums(pmat_gene_Adj<0.05)>0,] ), var='Module'),
						  geneCor_Pearson_Pval_AdjNetwork = tibble::rownames_to_column( as.data.frame( pmat_gene_Adj[rowSums(pmat_gene_Adj<0.05)>0,] ), var='Module'),

						  geneCor_Pearson_R_QsvaNetwork = tibble::rownames_to_column( as.data.frame( cormat_gene_qsva[rowSums(pmat_gene_qsva<0.05)>0,] ), var='Module'), 
						  geneCor_Pearson_Pval_QsvaNetwork = tibble::rownames_to_column( as.data.frame( pmat_gene_qsva[rowSums(pmat_gene_qsva<0.05)>0,] ), var='Module'),
						  
						  exonCor_Pearson_R_QsvaNetwork = tibble::rownames_to_column( as.data.frame( cormat_exon_qsva[rowSums(pmat_exon_qsva<0.05)>0,] ), var='Module') ,
						  exonCor_Pearson_Pval_QsvaNetwork = tibble::rownames_to_column( as.data.frame( pmat_exon_qsva[rowSums(pmat_exon_qsva<0.05)>0,] ), var='Module') ,
						  
						  jxnCor_Pearson_R_QsvaNetwork = tibble::rownames_to_column( as.data.frame( cormat_jxn_qsva[rowSums(pmat_jxn_qsva<0.05)>0,] ), var='Module') ,
						  jxnCor_Pearson_Pval_QsvaNetwork = tibble::rownames_to_column( as.data.frame( pmat_jxn_qsva[rowSums(pmat_jxn_qsva<0.05)>0,] ), var='Module') ),
						  file = '/dcl01/lieber/ajaffe/Steve/SNX19/tables/WGCNA_PolyA_DLPFC_Gene_Exon_Junction_Significant_EigengeneCorrelation_SNX19_Feature.xlsx')

################
# associations #
################

# gene set associations
library(clusterProfiler)
library(org.Hs.eg.db)

# split genes into modules, dropping grey
moduleGeneList_qsva = split(geneMap$EntrezID, netQsva$colors)
moduleGeneList_qsva = lapply(moduleGeneList_qsva, function(x) x[!is.na(x)])
moduleGeneList_qsva = moduleGeneList_qsva[-1]
moduleGeneList_adj = split(geneMap$EntrezID, netAdj$colors)
moduleGeneList_adj = lapply(moduleGeneList_adj, function(x) x[!is.na(x)])
moduleGeneList_adj = moduleGeneList_adj[-1]

## set universe of expressed genes
geneUniverse = as.character(geneMap$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

############################## 
## run enrichment analysis ###
##############################

## qSVA
goBP_qSVA <- compareCluster(moduleGeneList_qsva, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
goMF_qSVA <- compareCluster(moduleGeneList_qsva, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
goCC_qSVA <- compareCluster(moduleGeneList_qsva, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
kegg_qSVA <- compareCluster(moduleGeneList_qsva, fun = "enrichKEGG",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)

## and adjusted
goBP_Adj <- compareCluster(moduleGeneList_adj, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
goMF_Adj <- compareCluster(moduleGeneList_adj, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
goCC_Adj <- compareCluster(moduleGeneList_adj, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 1, qvalueCutoff  = 1,
				readable= TRUE)
