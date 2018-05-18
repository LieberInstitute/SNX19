library(tidyr)
library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size=14) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="bottom"	)) 
## SNX19 make plots
load('/dcl01/lieber/ajaffe/Steve/SNX19/eqtl_runs/Merge/rdas/FinalObject_Annotated_5LIBD-CMC_set_Merged_filtered_noLIBD_PGC_Missing_1e5sigFilt_reOrdered_morePGC.rda')

## load features of interest
#foi = read.csv('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_100517.csv',stringsAsFactors=FALSE)
#foi$Tag[foi$Tag=="junc8c.8d\xa0"] = "junc8c.8d"#Fixing Tag which is misread
#cisMerge$Tag = foi[match(cisMerge$Feature, foi$Feature), 'Tag']

load('/dcl01/lieber/ajaffe/Steve/SNX19/SNX19_8Features_of_Interest_crossBuild_03_22_18.rda')
cisMerge$Tag = foi[match(cisMerge$Feature, foi$Feature), 'Tag']
##
max(matrixStats ::rowMins(as.matrix(cisMerge[,grep("pvalue",colnames(cisMerge))]),na.rm=T))

#library(RColorBrewer)
#n <- 10
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))



## Reformat data for plotting
load('/dcl01/lieber/ajaffe/Steve/SNX19/Data/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_n495_updateMap_03082017.rda')
dat3 <- tidyr::gather(cisMerge, key=Dataset, value = Pvalue, DLPFC_PolyA_Linear_All_pvalue, DLPFC_RiboZero_Linear_All_pvalue, Caudate_Linear_All_pvalue, Hippo_Linear_All_pvalue, CMC_DLPFC_Linear_All_pvalue)

#Fixing dataset levels
dat3$Dataset <- factor(dat3$Dataset, levels=c('DLPFC_PolyA_Linear_All_pvalue','DLPFC_RiboZero_Linear_All_pvalue','CMC_DLPFC_Linear_All_pvalue','Hippo_Linear_All_pvalue',"Caudate_Linear_All_pvalue" ) )
#Fixing symbol stuff
dat3$Symbol[dat3$Symbol==""] <- NA
kk=match(unique(dat3$Symbol)[!is.na(unique(dat3$Symbol))], geneMap$Symbol)
dat3$Symbol <- factor(dat3$Symbol, geneMap$Symbol[kk][order(geneMap[kk,"Start"])]   )  
### 2MB region plot Alld
Region_features_all_snps_plot <- ggplot(data=subset(dat3, !is.na(dat3$Symbol)), aes(x=((exprsStart+exprsEnd)/2)/1e6, y= -log10(Pvalue),col=Symbol ) )
Region_features_all_snps_plot <- Region_features_all_snps_plot + facet_wrap(~Dataset,scales='free',nrow=1, labeller=labeller(
								 Dataset=c(DLPFC_PolyA_Linear_All_pvalue="DLPFC PolyA", 
										   DLPFC_RiboZero_Linear_All_pvalue="DLPFC Ribo", 
										   Caudate_Linear_All_pvalue="Caudate",
										   Hippo_Linear_All_pvalue = "Hippo",
										   CMC_DLPFC_Linear_All_pvalue = "CMC DLPFC" )) ) +
  geom_point() + 
  labs(x="Expressed Feature Chromosomal Position (Mb)", 
       y="-log10 (eQTL p value)") +   
  guides(col=guide_legend(title="Gene Symbol"))  +
  theme(legend.background = element_rect(colour = "black")) + 
  scale_y_continuous(limits = c(-log10(1e-4), 100),expand=c(0,0)) + theme_bw(base_size=14) +
  theme(legend.justification = c(1, 0.5), 
        legend.position = 'right', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(angle=45,size=10,hjust=1) ) 
ggsave(Region_features_all_snps_plot, filename='/dcl01/lieber/ajaffe/Steve/SNX19/plots/Region_features_all_snps_plot_allPvalues.pdf', width=11, height=8.5)

### 2MB region plot filtered
Region_features_all_snps_plot <- ggplot(data=subset(dat3, !is.na(dat3$Symbol) & dat3$pgc52_P <5e-8), aes(x=((exprsStart+exprsEnd)/2)/1e6, y= -log10(Pvalue),col=Symbol ) )
Region_features_all_snps_plot <- Region_features_all_snps_plot + facet_wrap(~Dataset,scales='free',nrow=1, labeller=labeller(
								 Dataset=c(DLPFC_PolyA_Linear_All_pvalue="DLPFC PolyA", 
										   DLPFC_RiboZero_Linear_All_pvalue="DLPFC Ribo", 
										   Caudate_Linear_All_pvalue="Caudate",
										   Hippo_Linear_All_pvalue = "Hippo",
										   CMC_DLPFC_Linear_All_pvalue = "CMC DLPFC" )) ) +
  geom_point() + 
  labs(x="Expressed Feature Chromosomal Position (Mb)", 
       y="-log10 (eQTL p value)") +   
  guides(col=guide_legend(title="Gene Symbol"))  +
  theme(legend.background = element_rect(colour = "black")) + 
  scale_y_continuous(limits = c(-log10(1e-4), 45),expand=c(0,0)) + theme_bw(base_size=14) +
  theme(legend.justification = c(1, 0.5), 
        legend.position = 'right', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x  = element_text(angle=45,size=10,hjust=1) ) 
ggsave(Region_features_all_snps_plot, filename='/dcl01/lieber/ajaffe/Steve/SNX19/plots/Region_features_all_snps_plot_all_PGCSig_Pvalues.pdf', width=11, height=8.5)

### SNX19 plot
dat3$Tag[is.na(dat3$Tag)] ="Other SNX19 Feature"
dat3$Tag <- factor(dat3$Tag, levels=c(foi$Tag ,"Other SNX19 Feature") )


SNX19_features_all_snps_plot <- ggplot(data=subset(dat3, Symbol=="SNX19"), aes(x=( (exprsStart+exprsEnd)/2)/1e6, y= -log10(Pvalue),col=Tag ) )
SNX19_features_all_snps_plot <- SNX19_features_all_snps_plot + facet_wrap(~Dataset,scales='free',nrow=1,labeller=labeller(
								 Dataset=c(DLPFC_PolyA_Linear_All_pvalue="DLPFC PolyA", 
										   DLPFC_RiboZero_Linear_All_pvalue="DLPFC Ribo", 
										   Caudate_Linear_All_pvalue="Caudate",
										   Hippo_Linear_All_pvalue = "Hippo",
										   CMC_DLPFC_Linear_All_pvalue = "CMC DLPFC" )) )+
  geom_point(shape=1,
             size=1.5,
             alpha=.7) + 
  labs(x="Expressed Feature Chromosomal Position (Mb)", 
       y="-log10 ( eQTL p value)") +
  guides(colour = guide_legend(override.aes = 
                                list(size = 3, shape=c(16),alpha=1 )))+
									 theme(legend.background = element_rect(colour = "black"))  + 
  theme(legend.background = element_rect(colour = "black")) + 
  theme(legend.justification = c(1, 0.5), 
        legend.position = 'right', 
        legend.title=element_blank(), #no legend title
        panel.grid.major = element_blank(), #no grid lines
        panel.grid.minor = element_blank(),
		axis.text.x  = element_text(angle=45,size=10,hjust=1)		) + 		#no grid lines
		scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") +
 # scale_colour_manual(values=c("blue","grey","green","purple","red","black","orange","teal")) +
  geom_point(data=subset(dat3, Symbol=="SNX19" &`pgc52_P`<5e-8), 
             col= "red",
             size=1,
             shape=4) + facet_wrap(~Dataset, scales='free',nrow=1,labeller=labeller(
								 Dataset=c(DLPFC_PolyA_Linear_All_pvalue="DLPFC PolyA", 
										   DLPFC_RiboZero_Linear_All_pvalue="DLPFC Ribo", 
										   Caudate_Linear_All_pvalue="Caudate",
										   Hippo_Linear_All_pvalue = "Hippo",
										   CMC_DLPFC_Linear_All_pvalue = "CMC DLPFC" )) ) +   scale_y_continuous(limits = c(0,100), expand = c(0, 0))

ggsave(SNX19_features_all_snps_plot, filename='/dcl01/lieber/ajaffe/Steve/SNX19/plots/SNX19_Feature_P_byDataset.pdf', width=11, height=8.5)

### SNX19 gwas positive plot
dat3$Tag[is.na(dat3$Tag)] ="Other SNX19 Feature"
#dat3$Tag = factor(dat3$Tag, levels = dat3$Tag[!duplicated(dat3$Tag)][order(rowSums(dat3[!duplicated(dat3$Tag),c('exprsStart','exprsEnd')])/2)] )

SNX19_features_all_snps_plot <- ggplot(data=subset(dat3, Symbol=="SNX19" & dat3$pgc52_P<5e-8 ), aes(x=( (exprsStart+exprsEnd)/2)/1e6, y= -log10(Pvalue),col=Tag ) )
SNX19_features_all_snps_plot <- SNX19_features_all_snps_plot + facet_wrap(~Dataset,scales='free',nrow=1,labeller=labeller(
								 Dataset=c(DLPFC_PolyA_Linear_All_pvalue="DLPFC PolyA", 
										   DLPFC_RiboZero_Linear_All_pvalue="DLPFC Ribo", 
										   Caudate_Linear_All_pvalue="Caudate",
										   Hippo_Linear_All_pvalue = "Hippo",
										   CMC_DLPFC_Linear_All_pvalue = "CMC DLPFC" )) )+
  geom_point(shape=1,
             size=1.5,
             alpha=1) + 
  labs(x="Expressed Feature Chromosomal Position (Mb)", 
       y="-log10 ( eQTL p value)") +
  guides(colour = guide_legend(override.aes = 
                                list(size = 3, shape=c(16),alpha=1 )))+
									 theme(legend.background = element_rect(colour = "black")) + 
  theme(legend.justification = c(1, 0.5), 
        legend.position = 'right', 
        legend.title=element_blank(), #no legend title
        panel.grid.major = element_blank(), #no grid lines
        panel.grid.minor = element_blank(),
		axis.text.x  = element_text(angle=45,size=10,hjust=1)		) + 		#no grid lines
		scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") + facet_wrap(~Dataset, scales='free',nrow=1,labeller=labeller(
								 Dataset=c(DLPFC_PolyA_Linear_All_pvalue="DLPFC PolyA", 
										   DLPFC_RiboZero_Linear_All_pvalue="DLPFC Ribo", 
										   Caudate_Linear_All_pvalue="Caudate",
										   Hippo_Linear_All_pvalue = "Hippo",
										   CMC_DLPFC_Linear_All_pvalue = "CMC DLPFC" )) ) +   scale_y_continuous(limits = c(0,45), expand = c(0, 0))  

ggsave(SNX19_features_all_snps_plot, filename='/dcl01/lieber/ajaffe/Steve/SNX19/plots/SNX19_Feature_P_byDataset_gwasPos.pdf', width=11, height=8.5)

### SNX19 scattergram
SNX19_8features_all_snps_plot <- ggplot(data=subset(dat3, Symbol=="SNX19" & Tag !="Other SNX19 Feature"), aes(x=Tag, y= -log10(Pvalue),col=Tag ) )
SNX19_8features_all_snps_plot <- SNX19_features_all_snps_plot + facet_wrap(~Dataset,scales='free',nrow=1,labeller=labeller(
								 Dataset=c(DLPFC_PolyA_Linear_All_pvalue="DLPFC PolyA", 
										   DLPFC_RiboZero_Linear_All_pvalue="DLPFC Ribo", 
										   Caudate_Linear_All_pvalue="Caudate",
										   Hippo_Linear_All_pvalue = "Hippo",
										   CMC_DLPFC_Linear_All_pvalue = "CMC DLPFC" )) )+
  geom_point(shape=1,
             size=1.5,
             alpha=1) + 
  labs(x="Expressed Feature Chromosomal Position (Mb)", 
       y="-log10 (eQTL p value)") +
  theme(legend.background = element_rect(colour = "black")) + 
  theme(legend.justification = c(1, 0.5), 
        legend.position = 'none', 
        legend.title=element_blank(), #no legend title
        panel.grid.major = element_blank(), #no grid lines
        panel.grid.minor = element_blank(),
		axis.text.x  = element_text(angle=45,size=10,hjust=1),
		axis.title.x=element_blank()		) + 	
		scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") +guides(colour = guide_legend(override.aes = 
                                list(size = 3, shape=c(16),alpha=1 )))+
  geom_point(data=subset(dat3, Symbol=="SNX19" &`pgc52_P`<5e-8 & Tag !="Other SNX19 Feature"), 
             col= "black",
             size=1,
             shape=4) + facet_wrap(~Dataset, scales='free',nrow=1,labeller=labeller(
								 Dataset=c(DLPFC_PolyA_Linear_All_pvalue="DLPFC PolyA", 
										   DLPFC_RiboZero_Linear_All_pvalue="DLPFC Ribo", 
										   Caudate_Linear_All_pvalue="Caudate",
										   Hippo_Linear_All_pvalue = "Hippo",
										   CMC_DLPFC_Linear_All_pvalue = "CMC DLPFC" )) ) +   scale_y_continuous(limits = c(0,100), expand = c(0, 0))

ggsave(SNX19_8features_all_snps_plot, filename='/dcl01/lieber/ajaffe/Steve/SNX19/plots/SNX19_8Feature_P_scattergram_byDataset.pdf', width=11, height=8.5)

### SNX19 scattergram gwasPos
dat3$Tag[is.na(dat3$Tag)] ="Other SNX19 Feature"

SNX19_8features_all_snps_plot <- ggplot(data=subset(dat3, Symbol=="SNX19" & Tag !="Other SNX19 Feature" & pgc52_P<5e-8), aes(x=Tag, y= -log10(Pvalue),col=Tag ) )
SNX19_8features_all_snps_plot <- SNX19_8features_all_snps_plot + facet_wrap(~Dataset,scales='free',nrow=1,labeller=labeller(
								 Dataset=c(DLPFC_PolyA_Linear_All_pvalue="DLPFC PolyA", 
										   DLPFC_RiboZero_Linear_All_pvalue="DLPFC Ribo", 
										   Caudate_Linear_All_pvalue="Caudate",
										   Hippo_Linear_All_pvalue = "Hippo",
										   CMC_DLPFC_Linear_All_pvalue = "CMC DLPFC" )) )+
  geom_point(shape=1,
             size=1.5,
             alpha=1) + 
  labs(x="Expressed Feature Chromosomal Position (Mb)", 
       y="-log10 (eQTL p value)") +
  theme(legend.background = element_rect(colour = "black")) + 
  theme(legend.justification = c(1, 0.5), 
        legend.position = 'none', 
        legend.title=element_blank(), #no legend title
        panel.grid.major = element_blank(), #no grid lines
        panel.grid.minor = element_blank(),
		axis.text.x  = element_text(angle=45,size=10,hjust=1),
		axis.title.x=element_blank()		) + 	
		scale_colour_brewer(palette = "Set3") +
		scale_fill_brewer(palette = "Set3") +   scale_y_continuous(limits = c(0,45), expand = c(0, 0))+guides(colour = guide_legend(override.aes = 
                                list(size = 3, shape=c(16),alpha=1 )))

ggsave(SNX19_8features_all_snps_plot, filename='/dcl01/lieber/ajaffe/Steve/SNX19/plots/SNX19_8Feature_P_scattergram_byDataset_gwasPos.pdf', width=11, height=8.5)


### Lifespan plots			 

#### DLPFC PolyA #### 
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_polyA/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_n495_updateMap_03082017.rda",verbose=T)
pd_DlpfcPolya = pd

#expression object
geneRpkm_DlpfcPolya = geneRpkm
geneMap_DlpfcPolya = geneMap

exonRpkm_DlpfcPolya = exonRpkm
exonMap_DlpfcPolya = exonMap

jRpkm_DlpfcPolya = jRpkm
jMap_DlpfcPolya = jMap

yExprs_DlpfcPolya = rbind(geneRpkm_DlpfcPolya, exonRpkm_DlpfcPolya, jRpkm_DlpfcPolya )
snpAll_DlpfcPolyA = snpAll
#### DLPFC RiboZero #### 
load("/users/ajaffe/Lieber/Projects/RNAseq/SNX19/DLPFC_RiboZero/rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_RiboZero_n485_updateMap.rda")
pd_DlpfcRibo = pd

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
#expression object
geneRpkm_CMC = geneRpkm
geneMap_CMC = geneMap

exonRpkm_CMC = exonRpkm
exonMap_CMC = exonMap

jRpkm_CMC = jRpkm
jMap_CMC = jMap

yExprs_CMC = rbind(geneRpkm_CMC, exonRpkm_CMC, jRpkm_CMC)
#######################			 
library(jaffelab)
#DLPFC PolyA
mod = model.matrix(~age + snpPC1 +snpPC2 + snpPC3 + RIN + totalAssignedGene +Sex, data=pd_DlpfcPolya) #forming the model matrix
cleanRpkm_DlpfcPolyA = t(jaffelab::cleaningY(y=log2(yExprs_DlpfcPolya[as.character(foi$Feature),]+1),mod=mod, P=2))
pdf("/dcl01/lieber/ajaffe/Steve/SNX19/plots/Ten_Feature_DLPFC_PolyA_agePlotter_lifespan.pdf", width=11,height=8.5)
 for (i in as.character(unique(foi$Feature)) ) {
jaffelab::agePlotter(cleanRpkm_DlpfcPolyA[,i],pd_DlpfcPolya$age,mainText=paste0(foi$Tag[match(i, foi$Feature)], "   LIBD DLPFC PolyA" ), smoothIt = FALSE)
	  }
dev.off()

#DLPFC RiboZero
mod = model.matrix(~age + snpPC1 +snpPC2 + snpPC3 + RIN + totalAssignedGene +Sex, data=pd_DlpfcRibo) #forming the model matrix
cleanRpkm_DlpfcRibo = t(jaffelab::cleaningY(y=log2(yExprs_DlpfcRibo[as.character(foi$Feature),]+1),mod=mod, P=2))
pdf("/dcl01/lieber/ajaffe/Steve/SNX19/plots/Ten_Feature_DLPFC_Ribo_agePlotter_lifespan.pdf", width=11,height=8.5)
 for (i in as.character(unique(foi$Feature)) ) {
jaffelab::agePlotter(cleanRpkm_DlpfcRibo[,i],pd_DlpfcRibo$age,mainText=paste0(foi$Tag[match(i, foi$Feature)], "   LIBD DLPFC RiboZero"  ), smoothIt = FALSE)
	  }
dev.off()

# Caudate #NO FETAL SAMPLES
#mod = model.matrix(~age + snpPC1 +snpPC2 + snpPC3 + RIN + totalAssignedGene +Sex, data=pd_Caudate) #forming the model matrix
#cleanRpkm_Caudate = t(jaffelab::cleaningY(y=log2(yExprs_Caudate[as.character(foi$Feature),]+1),mod=mod, P=2))
##pdf("/dcl01/lieber/ajaffe/Steve/SNX19/plots/Eight_Feature_Caudate_agePlotter_lifespan.pdf", width=11,height=8.5)
# for (i in as.character(unique(foi$Feature)) ) {
#jaffelab::agePlotter(cleanRpkm_Caudate[,i],pd_Caudate$age,mainText=paste0(foi$Tag[match(i, foi$Feature)] ), smoothIt = FALSE)
#	  }
#dev.off()

# Hippo
mod = model.matrix(~Age + snpPC1 +snpPC2 + snpPC3 + RIN + totalAssignedGene +Sex, data=pd_Hippo) #forming the model matrix
cleanRpkm_Hippo = t(jaffelab::cleaningY(y=log2(yExprs_Hippo[as.character(foi$Feature),]+1),mod=mod, P=2))
pdf("/dcl01/lieber/ajaffe/Steve/SNX19/plots/Ten_Feature_Hippo_agePlotter_lifespan.pdf", width=11,height=8.5)
 for (i in as.character(unique(foi$Feature)) ) {
jaffelab::agePlotter(cleanRpkm_Hippo[,i],pd_Hippo$Age ,mainText=paste0(foi$Tag[match(i, foi$Feature)], "   LIBD Hippocampus"  ), smoothIt = FALSE)
	  }
dev.off()

# get some numbers
length(unique(cisMerge[which(cisMerge$pgc52_P<5e-8),'SNP']))
length(unique(cisMerge[which(cisMerge$Symbol=="SNX19"),'SNP']))
length(unique(cisMerge[,'SNP']))


#DLPFC PolyA
mod = model.matrix(~age + snpPC1 +snpPC2 + snpPC3 + RIN + totalAssignedGene +Sex, data=pd_DlpfcPolya) #forming the model matrix
cleanRpkm_DlpfcPolyA = t(jaffelab::cleaningY(y=log2(yExprs_DlpfcPolya[as.character(foi$Feature),]+1),mod=mod, P=2))
dat = cbind(pd_DlpfcPolya, cleanRpkm_DlpfcPolyA )

dat = cbind(pd_DlpfcPolya, t(log2(yExprs_DlpfcPolya[as.character(foi$Feature),]+1)) )

colnames(dat)[42:51] <- foi[ match(colnames(dat)[42:51], foi$Feature ), 'Tag']
dat = cbind(dat, rs948085 = t(snpAll_DlpfcPolyA["rs948085",] ) )
dat$rs948085 = as.factor(dat$rs948085)
dat$Race_genotype <- paste(dat$Race, dat$rs948085, sep="_")

library(ggplot2)
ggplot(data=dat, aes_string(x='age', y='junc8.10', col='rs948085' ) ) + geom_point()

rs10791097:130718630:T:G
