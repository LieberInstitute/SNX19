#qsub -l bluejay -l mf=15G,h_vmem=18G,h_stack=256M -pe local 9 -R y -M stephen.semick@gmail.com -cwd -b y R CMD BATCH SNX19_GTEX_Analysis.R
###
library(GenomicRanges)
library(matrixStats)
library(MatrixEQTL)
library(jaffelab)
library(parallel)
##########################
### load GTEx data
load('/users/ajaffe/Lieber/Projects/RNAseq/SNX19/GTEX/rawAndRpkmCounts_plusGenotype_SNX19_GTEX_n8206.rda')

#load annotation information
load("/dcl01/lieber/ajaffe/Steve/SNX19/Data/ensembl_V75_feature_to_Tx.rda") #load in transcript data (object is "allTx")


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


### get indices to split by region ####
rIndexes=split(1:nrow(pd), pd$SMTSD)
lengths(rIndexes)
#names(rIndexes) = gsub(" ", "_", ss(ss(names(rIndexes), 
#	"Brain - ", 2), " \\("))
	
########################### 
## analysis by region #####
###########################

inds = seq(along=rIndexes)
names(inds) = names(rIndexes)
inds=inds[1:2]
## gene level
SNX19_allFeature_EqtlGtex = mclapply(inds, function(i) {
	ii = rIndexes[[i]]
	### Tissue specific PCA
	pcaexprs3features = prcomp(t(exprs3features[,ii]+1)) #principle component decomposition on the transpose of that matrix + 1
	all3PCs = pcaexprs3features$x[,1:15] #take the first 15 PCs from that matrix
	
	
	## model
	mod = model.matrix(~pd[ii,'snpPC1'] + pd[ii,'snpPC2']  + pd[ii,'snpPC3'] + pd[ii,'snpPC4'] + pd[ii,'snpPC5'] + all3PCs) #forming the model matrix from the SNP PCs and the expression PCs
	colnames(mod)[2:ncol(mod)] = c(paste0("snpPC",1:5), paste0("exprsPC",1:15)) #Renaming the column names of the model matrix

	covs = SlicedData$new(t(mod[,-1])) #This part employs the "MatrixEQTL" package
	
	##SNP formatting
	snpSlice = SlicedData$new(snp2[,ii]) #formating the snp data for the MatrixEQTL package 
	snpSlice$ResliceCombined(sliceSize = 5000) 	

	## SNP Position
	snpspos = snpMapAll[,c("SNP","CHR","POS")] #snpspos is just a subset of snpMapAll containing the variables SNP, CHR, and POS
	snpspos$CHR = paste0("chr",snpspos$CHR) #concatenating the string "chr" with the numbers in "CHR"
	colnames(snpspos) = c("name","chr","pos") #renaming the snppos columns

	## Position objects
	posGene = geneMap[,c('Chr', 'Start', 'End')] #gene position information taken from the geneMap
	posGene$name = rownames(geneMap) #gene position names variable made from rownames of geneMap
	posGene = posGene[,c('name','Chr', 'Start', 'End')] #reordering the columns of the posGene object

	posExon = exonMap[,c('Chr', 'Start', 'End')] #exon position information taken from exonMap
	posExon$name = rownames(exonMap) #exon position names variable made from rownames of exonMap
	posExon = posExon[,c('name', 'Chr', 'Start', 'End')] #reordering the columns of the posExon dataframe

	posJxn = as.data.frame(jMap)[,c('seqnames','start','end')] #extracting junction position information from the jMap object after converting it to a dataframe
	posJxn$name = names(jMap) #junction position names are created from the names of the jMap object
	posJxn = posJxn[,c('name','seqnames','start','end')] #reordering the columns of posJxn
	posJxn <- plyr::rename(posJxn,c("seqnames"="Chr", 
                         'start' = 'Start', 
                           'end'='End')) #rename "chrpos" to "snp_coord"
						   
	### Combining objects

	exprsSlice = SlicedData$new(exprs3features[,ii])
	exprsSlice$ResliceCombined(sliceSize = 5000)
	pos = rbind(posGene, posExon, posJxn) #combining position information for the three different types of expression data


	meJoint = Matrix_eQTL_main(snps=snpSlice,          #SNP information
                           gene = exprsSlice,      # Gene information
                           cvrt = covs,            #Covariate information
                           output_file_name.cis = paste0("/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/scratch/",names(rIndexes)[i], ".txt"),
                           pvOutputThreshold.cis = 1,  
                           pvOutputThreshold=0,
                           snpspos = snpspos, 
                           genepos = pos, 
                           useModel = modelLINEAR,  # modelANOVA or modelLINEAR or modelLINEAR_CROSS
                           pvalue.hist = 100, min.pv.by.genesnp = TRUE, cisDist=1e6)
##### Eqtl object						   
eqtl = meJoint$cis$eqtls

##### Annotating						   
colnames(eqtl)[2] = "feature" #adding a column called feature
eqtl$feature = as.character(eqtl$feature) #the feature column contains expression information
eqtl$snps = as.character(eqtl$snps) #converting the eqtl snps variable to a character variable
colnames(eqtl)[colnames(eqtl) %in% c('statistic','pvalue','FDR','beta')] = paste0("All_", colnames(eqtl)[colnames(eqtl) %in% c('statistic','pvalue','FDR','beta')]) #Adding "ALL_ to "statistics, p-value, FDR, and beta" to help distinguish from future analyses 
eqtl$UniqueID = paste0(eqtl$snps, ".", eqtl$feature) #making a variable called Unique ID... Just a combination of SNP and feature
eqtl = eqtl[,c('UniqueID', 'snps', 'feature', 'All_statistic', 'All_pvalue', 'All_FDR', 'All_beta')] #just reordering the columns of the dataframe
	return(eqtl)
}, mc.cores=8, mc.preschedule=F)


save(SNX19_allFeature_EqtlGtex, file="/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/rdas/SNX19_allFeature_EqtlGtex.rda", compress=TRUE) #saving dlpfc object in a .rda file					

