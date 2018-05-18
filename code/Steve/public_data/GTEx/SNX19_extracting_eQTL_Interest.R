load("/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/rdas/SNX19_allFeature_EqtlGtex.rda")
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



##################### Proccess eQTL Results


#Throw out failed datasets
SNX19_allFeature_EqtlGtex <- SNX19_allFeature_EqtlGtex[sapply(SNX19_allFeature_EqtlGtex,is.data.frame)]

#define matching reference
ref_UniqueID=SNX19_allFeature_EqtlGtex[[1]]$UniqueID

#
SNX19_allFeature_EqtlGtex_pca_prep <- lapply (names(SNX19_allFeature_EqtlGtex), function(x) { 

#Appending name of study to column names
y = SNX19_allFeature_EqtlGtex[[x]] 
colnames(y)[colnames(y)%in%c("All_statistic", "All_pvalue", "All_FDR", "All_beta")] <- paste0(x,"_",colnames(y)[colnames(y)%in%c("All_statistic", "All_pvalue", "All_FDR", "All_beta")] )

#Lining up all of eQTL Results 
y = y[match(ref_UniqueID,y$UniqueID), ]
return(y)
} )
names(SNX19_allFeature_EqtlGtex_pca_prep) <-  names(SNX19_allFeature_EqtlGtex)


###
reoriented_eqtl_res <- lapply(SNX19_allFeature_EqtlGtex_pca_prep[2:length(SNX19_allFeature_EqtlGtex_pca_prep)], function(x) {
x[,!colnames(x)%in%c("UniqueID", "snps", "feature", "All_beta")]
})

#Extracting SNPs of interest from GTEX eQTL
bound <- do.call("cbind", reoriented_eqtl_res)
bound <- cbind(SNX19_allFeature_EqtlGtex_pca_prep[[1]] , bound)
colnames(bound) <- gsub(".*\\.","",colnames(bound))

snp_of_interest=c('rs10791097',
'rs34529622',
'rs3794133',
'rs3831404',
'rs4337054',
'rs73028899',
'rs78968418',
'rs10791114')

### Intermediate Saving
save(bound,file='/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/rdas/gtex_results_wide_raw.rda')
### Annotating
bound$eqtlType = "Jxn" #fill the variable eqtlType with "Jxn"
bound$eqtlType[bound$feature %in% rownames(exonMap)] = "Exon" #when the feature variable starts with e call this an exon
bound$eqtlType[grep("^ENS", bound$feature)] = "Gene" #when the feature variable starts with ENS call this a gene
#
###### SNP Annotation
snpAnno = snpMapAll[match(bound$snps, snpMapAll$SNP), names(snpMapAll)[!names(snpMapAll) %in% "X.C.M"] ] #extracting annotation information for each SNP in the eqtl$snps
bound$snpRsNum = gsub("\\:.*","",bound$snpRsNum)

##Dropping the X.C.M variable in snpMapAll
snpAnno$CHR = paste0("chr", snpAnno$CHR) 
snpAnno <- plyr::rename(snpAnno, c("CHR"="snpChr", #renaming variables in the snpAnno object
                  "POS"="snpPos", 
                  "COUNTED"="snpCountedAllele",
                  "ALT"="snpRefAllele",
                  "chrpos" = "snpChrPos",
                  "SNP" ="snpRsNum"
                  ))
bound = cbind(bound, snpAnno) #adding the SNP annotation information into the dlpfc object
#
#### Feature Annotation
geneMap$Class = "InEns"
geneMap$Geneid = rownames(geneMap) #the gene id 
geneMap$Type <- "Gene"
exonMap$Class = "InEns"
exonMap$Type <- "Exon"
tmp = as.data.frame(jMap) #creating tmp by converting jMap to dataframe
tmp <- plyr::rename(tmp,  #renaming some of the variable names in the posJxn object
                       c("seqnames"="Chr", 
                         'start' = 'Start', 
                         'end'='End')) #rename "chrpos" to "snp_coord"
jMapDf = tmp[, c("Chr", "Start", "End", "code", 'newGeneID', "newGeneSymbol")] #creating jMapDf from some of the columns in tmp
jMapDf <- plyr::rename(jMapDf,  #renaming some of the variable names in the jMapDf object
                    c("code"="Class", 
                      'newGeneID' = 'Geneid', 
                      'newGeneSymbol'='Symbol')) #rename "chrpos" to "snp_coord"
jMapDf$Type <- "Jxn"
#
#
map = rbind(geneMap[,colnames(jMapDf)], exonMap[,colnames(jMapDf)], jMapDf) #Combining geneMap, exonMap, and jMap into a single Map
map$featureID = rownames(map)
#### Adding transcript information
##load tx info earlier
mmTx = match(map$featureID, names(allTx)) #matching transcripts with data in map. Match over ENSEMBL IDs
tx = CharacterList(vector("list", nrow(map))) #initializing an empty CharacterList
tx[!is.na(mmTx)] = allTx[mmTx[!is.na(mmTx)]] #using !is.na() to maintain correct indexing while filling the tx object with
map$NumTx = elementNROWS(tx)
map$WhichTx = sapply(tx, paste0, collapse=";")
#
#
####
colnames(map)[colnames(map)%in%c("Chr","Start","End","Class","Type","Geneid","Symbol","NumTx","WhichTx")] = paste0("exprs", colnames(map)[colnames(map)%in%c("Chr","Start","End","Class","Type","Geneid","Symbol","NumTx","WhichTx")] ) #adding "exprs" to the names of several columns
map = map[match(bound$feature, map$featureID),] #keeping only the rows of map that have corresponding features in eqtl
bound = cbind(bound, map[,!colnames(map) %in% 'featureID']) #adding all the variables in map to eqtl (except for featureID)
###### Distance metrics
dToStart = bound$snpPos - bound$exprsStart #distance from single nucleotide polymorphism to expression start site
dToEnd = bound$snpPos - bound$exprsEnd #distance from single nucleotide polymorphism to expression end site
bound$distSnpToFeature = apply(cbind(dToStart,dToEnd), 1, function(x) x[which.min(abs(x))]) #minimum distance
bound$snpRsNum <- gsub(':.*', '', bound$snpRsNum)

save(bound, file='/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/rdas/gtex_results_wide_annotated.rda')
SOI <- bound[bound$snpRsNum%in%snp_of_interest,]
write.csv(SOI,file='/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/csvs/snps_of_interest_gtex_subset_080417_unflipped_annotated.csv',row.names=FALSE)

boundSubset <- bound[(rowSums(bound[,grepl("pvalue",colnames(bound))]<1e-4)>0),]
write.csv(boundSubset, file='/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/csvs/sig_eqtl_1e4_gtex_subset_080417_unflipped_annotated.csv',row.names=FALSE)
save(boundSubset, file='/dcl01/lieber/ajaffe/Steve/SNX19/GTEX/rdas/sig_eqtl_1e4_gtex_subset_080417_unflipped_annotated.rda')

hist( rowSums(bound[,grepl("pvalue",colnames(bound))]<1e-4) )

