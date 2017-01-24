##############
## preprocess SNX19 Brain Counts
##############

### load scripts
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

junctionCount = function(junctionFiles, 
	sampleNames=names(junctionFiles), 
	strandSpecific = FALSE, illuminaStranded=FALSE,
	minCount = 1, maxCores=NULL) {
	
	require(GenomicRanges,quietly=TRUE)
	require(parallel,quietly=TRUE)

	if(is.null(maxCores)) {
		maxCores=1
	}
		
	names(junctionFiles) = sampleNames
	cat("Reading in data.\n")
	if(all(is.character(junctionFiles))) {
		theData = mclapply(junctionFiles, function(x) {
			y = read.delim(x, skip = 1, header=FALSE, 
				col.names = c("chr", "start","end", "strand", "count"), 
				colClasses = c("character", "integer", "integer", "character","integer"))
			y = y[y$count >= minCount,] # filter based on min number
			gr = GRanges(y$chr, IRanges(y$start, y$end), 
				strand=y$strand,count = y$count)
			return(gr)
		}, mc.cores=maxCores)
	} else {
		theData = junctionFiles
		stopifnot(all(sapply(theData, class)=="GRanges"))
	}
	cat("Creating master table of junctions.\n")

	## turn into GRangesList
	### THIS STEP IS SLOW...
	grList = GRangesList(theData)

	# each dataset should be checked
	if(illuminaStranded & strandSpecific) {
		grList = GRangesList(mclapply(grList, function(x) {
			strand(x) = ifelse(strand(x)=="+", "-","+")
			return(x)
		},mc.cores=maxCores))
	}
	
	## get into GRanges object of unique junctions
	fullGR = unlist(grList)
	if(!strandSpecific) strand(fullGR) = "*"
	
	fullGR = fullGR[!duplicated(fullGR)] # or unique(fullGR)
	fullGR = sort(fullGR)
	fullGR$count = NULL

	cat(paste("There are", length(fullGR), "total junctions.\n"))
	
	cat("Populating count matrix.\n")

	jNames = paste0(as.character(seqnames(fullGR)),":",
			start(fullGR),"-",end(fullGR),"(",as.character(strand(fullGR)),")")

	## match GRanges
	options(warn=-1)
	mList = mclapply(grList, match, fullGR, 
		ignore.strand = !strandSpecific, mc.cores=maxCores)
	options(warn=0)
	
	countList = mList # initiate 
	M = length(jNames)

	## fill in matrix
	for(i in seq(along=grList)) {
		if(i %% 25 == 0) cat(".")
		cc = rep(0,M)
		cc[mList[[i]]] = theData[[i]]$count
		countList[[i]] = Rle(cc)
	}
	countDF = DataFrame(countList, row.names=jNames)
	
	names(fullGR) = jNames
	## return matrix and GRanges object
	out = list(countDF = countDF, anno = fullGR)
	return(out)
}

### load packages
library(GenomicRanges)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(biomaRt)
library(parallel)
library(GenomicAlignments)

# load phenotype data
load("phenotype_annotated_SNX19_plusGenotype_LIBD_DLPFC.rda")

#########################
#### GENE COUNTS ########

### gene count files
geneFn = paste0("/dcl01/lieber/ajaffe/Brain/SNX19/Counts/Gene/",
	pd$RNum, "_SNX19_Ensembl75_genes_DLPFC.counts")
names(geneFn) = pd$RNum

### read in annotation ##
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]
rownames(geneMap) = geneMap$Geneid
geneMap$Chr = ss(geneMap$Chr, ";")
geneMap$Start = as.numeric(ss(geneMap$Start, ";"))
tmp = strsplit(geneMap$End, ";")
geneMap$End = as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand = ss(geneMap$Strand, ";")
geneMap$Geneid = NULL

### biomart 
ensembl = useMart("ENSEMBL_MART_ENSEMBL", # VERSION 75, hg19
	dataset="hsapiens_gene_ensembl",
	host="feb2014.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
	values=rownames(geneMap), mart=ensembl)
geneMap$Symbol = sym$hgnc_symbol[match(rownames(geneMap), sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene[match(rownames(geneMap), sym$ensembl_gene_id)]
	
## counts
geneCountList = mclapply(geneFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=12)
geneCounts = do.call("cbind", geneCountList)
rownames(geneCounts) = rownames(geneMap)
geneCounts = geneCounts[,pd$RNum] # put in order

# make RPKM
bg = matrix(rep(pd$totalMapped), nc = nrow(pd), 
	nr = nrow(geneCounts),	byrow=TRUE)
widG = matrix(rep(geneMap$Length), nr = nrow(geneCounts), 
	nc = nrow(pd),	byrow=FALSE)
geneRpkm = geneCounts/(widG/1000)/(bg/1e6)

# number of reads assigned
geneStatList = lapply(paste0(geneFn, ".summary"), 
	read.delim,row.names=1)
geneStats = do.call("cbind", geneStatList)
colnames(geneStats) = pd$RNum
pd$totalAssignedGene = as.numeric(geneStats[1,] / colSums(geneStats))

#################################
######### EXON COUNTS ###########

### exon count files
exonFn = paste0("/dcl01/lieber/ajaffe/Brain/SNX19/Counts/Exon/",
	pd$RNum, "_SNX19_Ensembl75_exons_DLPFC.counts")
names(exonFn) = pd$RNum

### read in annotation ##
exonMap = read.delim(exonFn[1], skip=1, as.is=TRUE)[,1:6]
exonMap$Symbol = sym$hgnc_symbol[match(exonMap$Geneid, sym$ensembl_gene_id)]
exonMap$EntrezID = sym$entrezgene[match(exonMap$Geneid, sym$ensembl_gene_id)]

## counts
exonCountList = mclapply(exonFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=12)
exonCounts = do.call("cbind", exonCountList)
exonCounts = exonCounts[,pd$RNum] # put in order

## remove duplicated
eMap = GRanges(exonMap$Chr, IRanges(exonMap$Start, exonMap$End))
keepIndex= which(!duplicated(eMap))
exonCounts = exonCounts[keepIndex,]
exonMap = exonMap[keepIndex,]

rownames(exonMap) = with(exonMap, paste0(Chr, 
	":",Start,"-", End, "(", Strand, ")"))
rownames(exonCounts) = rownames(exonMap)

## make RPKM
bgE = matrix(rep(pd$totalMapped), nc = nrow(pd), 
	nr = nrow(exonCounts),	byrow=TRUE)
widE = matrix(rep(exonMap$Length), nr = nrow(exonCounts), 
	nc = nrow(pd),	byrow=FALSE)
exonRpkm = exonCounts/(widE/1000)/(bgE/1e6)

########################
##### JUNCTIONS #############

countJunctionsFromBam <- function(bamFile) {
	cat(".")
	GA = readGAlignments(bamFile)
	GA_gapped = GA[njunc(GA) > 0]
	junc = unlist(junctions(GA_gapped))
	uniqueJunc = unique(junc)
	strand(uniqueJunc)="*"
	uniqueJunc = unique(uniqueJunc)
	counts = countOverlaps(uniqueJunc, junc, type="equal")
	uniqueJunc$count = counts
	return(uniqueJunc)
}

bamFiles = pd$bamFile
names(bamFiles)=pd$RNum
junctionCounts = mclapply(bamFiles, countJunctionsFromBam, mc.cores=8)

### get junction counts
juncCounts = junctionCount(junctionCounts,maxCores=12)

## annotate junctions
load("/users/ajaffe/Lieber/Projects/RNAseq/ensembl_hg19_v75_junction_annotation.rda")

anno = juncCounts$anno
seqlevels(anno, force=TRUE) = paste0("chr", c(1:22,"X","Y","M"))

## add additional annotation
anno$inEnsembl = countOverlaps(anno, theJunctions, type="equal") > 0
anno$inEnsemblStart = countOverlaps(anno, theJunctions, type="start") > 0
anno$inEnsemblEnd = countOverlaps(anno, theJunctions, type="end") > 0

oo = findOverlaps(anno, theJunctions, type="equal")
anno$ensemblGeneID = NA
anno$ensemblGeneID[queryHits(oo)] = as.character(theJunctions$ensemblID[subjectHits(oo)])
anno$ensemblSymbol = NA
anno$ensemblSymbol[queryHits(oo)] = theJunctions$symbol[subjectHits(oo)]
anno$ensemblStrand = NA
anno$ensemblStrand[queryHits(oo)] = as.character(strand(theJunctions)[subjectHits(oo)])
anno$ensemblTx = CharacterList(vector("list", length(anno)))
anno$ensemblTx[queryHits(oo)] = theJunctions$tx[subjectHits(oo)]
anno$numTx = elementNROWS(anno$ensemblTx)

## clean up gene symbol
anno$ensemblSymbol = geneMap$Symbol[match(anno$ensemblGeneID, rownames(geneMap))]

# sequence of acceptor/donor sites
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
left = right = anno
end(left) = start(left) +1
start(right) = end(right) -1

anno$leftSeq  = getSeq(Hsapiens, left)
anno$rightSeq = getSeq(Hsapiens, right)

## junction code
anno$code = ifelse(anno$inEnsembl, "InEns", 
	ifelse(anno$inEnsemblStart & anno$inEnsemblEnd, "ExonSkip",
	ifelse(anno$inEnsemblStart | anno$inEnsemblEnd, "AltStartEnd", "Novel")))

## b/w exons and junctions
exonGR = GRanges(exonMap$Chr,IRanges(exonMap$Start, exonMap$End))
anno$startExon = match(paste0(seqnames(anno),":",start(anno)-1), 
	paste0(seqnames(exonGR), ":", end(exonGR)))
anno$endExon = match(paste0(seqnames(anno),":",end(anno)+1),
	paste0(seqnames(exonGR), ":", start(exonGR)))
g = data.frame(leftGene = exonMap$Geneid[anno$startExon],
	rightGene = exonMap$Geneid[anno$endExon],
	leftGeneSym = exonMap$Symbol[anno$startExon],
	rightGeneSym = exonMap$Symbol[anno$endExon],
	stringsAsFactors=FALSE)
g$newGene = NA
g$newGeneSym = NA
g$newGene[which(g$leftGene==g$rightGene)] = 
	g$leftGene[which(g$leftGene==g$rightGene)] 
g$newGeneSym[which(g$leftGene==g$rightGene)] = 
	g$leftGeneSym[which(g$leftGene==g$rightGene)] 
g$newGene[which(g$leftGene!=g$rightGene)] = 
	paste0(g$leftGene,"-",g$rightGene)[which(g$leftGene!=g$rightGene)] 
g$newGeneSym[which(g$leftGene!=g$rightGene)] = 
	paste0(g$leftGeneSym,"-",g$rightGeneSym)[which(g$leftGene!=g$rightGene)] 
g$newGene[which(is.na(g$newGene) & is.na(g$leftGene))] = 
	g$rightGene[which(is.na(g$newGene) & is.na(g$leftGene))] 
g$newGene[which(is.na(g$newGene) & is.na(g$rightGene))] = 
	g$leftGene[which(is.na(g$newGene) & is.na(g$rightGene))] 
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] = 
	g$rightGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] 
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] = 
	g$leftGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] 
g$newGeneSym[g$newGeneSym==""] = NA
g$newGeneSym[g$newGeneSym=="-"] = NA
anno$newGeneID = g$newGene
anno$newGeneSymbol = g$newGeneSym
anno$isFusion = grepl("-", anno$newGeneID)

## extract out
jMap = anno
jCounts = juncCounts$countDF
jCounts = as.data.frame(jCounts[names(jMap),pd$RNum])

# MAPPED PER MILLION
mappedPerM = pd$totalMapped/1e6
countsM = DataFrame(mapply(function(x,d) x/d, jCounts , mappedPerM))
rownames(jCounts) = rownames(countsM) = names(jMap)
jRpkm = as.data.frame(countsM)

# fix SNPs that are in/dels 
library(jaffelab)
library(stringr)
ncRef= nchar(snpMapAll$ALT)
ncCount= nchar(snpMapAll$COUNTED)
snpMapAll$Type = "SNV"
snpMapAll$Type[ncRef > ncCount] = "Deletion"
snpMapAll$Type[ncRef < ncCount] = "Insertion"

snpMapAll$newRef = snpMapAll$ALT
snpMapAll$newCount = snpMapAll$COUNT

# deletion
dIndex = which(snpMapAll$Type=="Deletion")
snpMapAll$newRef[dIndex] = substr(snpMapAll$ALT[dIndex], 
	2,nchar(snpMapAll$ALT[dIndex]))
snpMapAll$newCount[dIndex] = "-"

# insertion
iIndex = which(snpMapAll$Type=="Insertion")
snpMapAll$newRef[iIndex] = "-"
snpMapAll$newCount[iIndex] = substr(snpMapAll$COUNTED[iIndex], 
	2,nchar(snpMapAll$COUNTED[iIndex]))
head(snpMapAll[snpMapAll$Type != "SNV",],10)

### classify other map
snpInfo$LIBD_newRef = snpMapAll$newRef[match(snpInfo$LIBD_SNP_ID,snpMapAll$SNP)]
snpInfo$LIBD_newCount = snpMapAll$newCount[match(snpInfo$LIBD_SNP_ID,snpMapAll$SNP)]

### save counts
save(pd, jMap, jCounts, geneCounts, geneMap, exonCounts, exonMap, 
	geneRpkm, exonRpkm, jRpkm,snpInfo, snpMat, 
	snpAll, snpMapAll, compress=TRUE,
	file="rawAndRpkmCounts_plusGenotype_SNX19_DLPFC_n495_updateMap.rda")
