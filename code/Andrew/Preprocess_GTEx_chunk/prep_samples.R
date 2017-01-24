#####
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)
library(stringr)

# phenotype data

# load data and filter
load("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/gtexPd.Rdata")
gtexPd$sra_accession = str_trim(gtexPd$sra_accession)
gtexPd = gtexPd[!is.na(gtexPd$sra_accession),]

## drop those w/o SRR
gtexPd = gtexPd[which(gtexPd$SAMPLE_USE == 
	"Seq_RNA_WTSS; Seq_RNA_Expression"),]
pd = gtexPd
pd$IndividualID = ss(pd$SUBJID, "-", 2)

############################
#### genotype data #########
bfile = "/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Genotypes/GTEX_5M_Genotype_Data_Imputed_maf005_geno10_hwe1e6"

fam = read.table(paste0(bfile, ".fam"),as.is=TRUE)
colnames(fam) = c("SAMPID", "SAMPID", "MID","PID","Sex","Pheno")
fam$IndividualID = ss(fam$SAMPID, "-", 2)
pd = pd[which(pd$IndividualID %in% fam$IndividualID),] # keep samples w/ genotypes

# fam = fam[!duplicated(fam$IndividualID),] # remove 650 if they have 1M
# famOut = fam[which(fam$IndividualID %in% pd$IndividualID),]
# write.table(famOut[,1:2], "samples_to_extract.txt",
	# col.names=FALSE, row.names=FALSE, quote=FALSE)

# ##### get common SNPs for MDS components
newbfile = "/dcl01/lieber/ajaffe/Brain/SNX19/GTEX/Plink/samples_with_RNAseq_GTEX"
# system(paste("plink --bfile", bfile, 
	# "--keep samples_to_extract.txt --make-bed --geno 0.1 --maf 0.05",
	# "--hwe 0.000001 --out", newbfile,"--memory 225000"))
# system(paste("plink --bfile", newbfile,"--indep 100 10 1.25 --out", newbfile))
# system(paste0("plink --bfile ", newbfile, 
	# " --cluster --mds-plot 5 --extract ",
	# newbfile, ".prune.in --out ", newbfile))

# write out RNums
# cat(pd$sra_accession, file="GTEX_Samples_SNX19_analysis.txt",sep="\n")

####################
### cut and quantify using shell script
####################

# add bam info
pd$bamFile = paste0("/dcl01/lieber/ajaffe/Brain/SNX19/GTEX/BAM/",
	pd$sra_accession, "_SNX19.bam")
all(file.exists(pd$bamFile)) #  TRUE

getTotalMapped = function(bamFile, mc.cores=1, returnM = TRUE) {
	thecall = paste("samtools idxstats",bamFile)
	tmp = parallel::mclapply(thecall, function(x) {
		cat(".")
		xx = system(x,intern=TRUE)
		xx = do.call("rbind", strsplit(xx, "\t"))
		d = data.frame(chr=xx[,1], L=xx[,2], mapped = xx[,3],
			stringsAsFactors=FALSE)
		d
	},mc.cores=mc.cores)
	
	out = list(totalMapped = sapply(tmp, function(x) sum(as.numeric(x$mapped[x$chr %in% paste0("chr", c(1:22,"X","Y"))]))),
		mitoMapped = sapply(tmp, function(x) as.numeric(x$mapped[x$chr=="chrM"])))
	return(out)
}

totMapped = getTotalMapped(pd$bamFile,mc.cores=12)
pd$totalMapped = totMapped$totalMapped

## drop two samples
pd = pd[pd$totalMapped > 500,]

######################
### get SNPs #########
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

pgc = read.delim("/users/ajaffe/Lieber/Projects/RNAseq/controlPaper/pgc/daner_PGC_SCZ52_0513a.gz.p4.clump.areator.sorted.1mhc.txt",
	as.is=TRUE)[18,] # SNX19 locus
tmp = unlist(strsplit(pgc$LD.friends.0.6..p0.001,","))
tmp2 = strsplit(pgc$LD.friends.0.6..p0.001,",")
theSnps = data.frame(name = c(pgc$SNP, ss(tmp, "\\(")),
	R2 = as.numeric(c(rep(1,nrow(pgc)), ss(ss(tmp, "\\(",2),"/"))),
	Dist = as.numeric(c(rep(0,nrow(pgc)), 
		gsub(")", "", ss(ss(tmp, "\\(",2),"/",2),fixed=TRUE))))
theSnps$name = as.character(theSnps$name)	

### get single SNP stats
awkCall = "awk '$1 == \"chr11\" {print $0}' /dcs01/ajaffe/Brain/SNPs/scz2.snp.results.txt"
headerCall = "awk NR==1 /dcs01/ajaffe/Brain/SNPs/scz2.snp.results.txt"
pgcStats = read.table(pipe(awkCall), as.is=TRUE,header=TRUE)
colnames(pgcStats) = read.table(pipe(headerCall),as.is=TRUE)[1,]

## match up
mm = match(theSnps$name, pgcStats$snpid)
theSnps = cbind(theSnps, pgcStats[mm,])
theSnps$chrpos = paste0(theSnps$hg19chrc, ":", theSnps$bp)

# ##### get plink name in LIBD
# bimFile = read.table(paste0(newbfile,".bim"), as.is=TRUE, 
	# col.names=c("CHR","NAME","CM","POS","COUNTED","REF"))
# bimFile$CHR = paste0("chr", bimFile$CHR)
# bimFile$chrpos = paste0(bimFile$CHR, ":", bimFile$POS)
# theSnps$GTEX_SNP_ID = bimFile$NAME[ # match
	# match(theSnps$chrpos,bimFile$chrpos)]

# cat(theSnps$GTEX_SNP_ID, file="snps_to_extract.txt",sep="\n")
# system(paste("plink --bfile", newbfile, 
	# "--keep samples_to_extract.txt --extract snps_to_extract.txt",
	# "--recode A-transpose --out snps_at_SNX19_GTEX --memory 225000"))

## read in genotypes
genotypes  = read.delim("snps_at_SNX19_GTEX.traw",as.is=TRUE)
snpMap = genotypes[,1:6]
snpMap$chrpos = paste0("chr", snpMap$CHR, ":", snpMap$POS)

genotypes = genotypes[,-(1:6)]
rownames(genotypes) = snpMap$SNP
colnames(genotypes) = ss(colnames(genotypes), "\\.",2)
snp = genotypes[,pd$IndividualID]

## put in order
mm2 = match(theSnps$chrpos, snpMap$chrpos)
theSnps = theSnps[!is.na(mm2),]
snpMap = snpMap[mm2[!is.na(mm2)],]
snp = snp[mm2[!is.na(mm2)],]
theSnps$GTEX_COUNTED = snpMap$COUNTED
theSnps$GTEX_REF = snpMap$ALT
rownames(snp) = rownames(theSnps) = theSnps$name

snpInfo = theSnps
snpMat = as.matrix(snp)

### add MDS
mds = read.table(paste0(newbfile, ".mds"), 
	header=TRUE,as.is=TRUE,row.names=1)[,-(1:2)]
colnames(mds) = paste0("snpPC",1:ncol(mds))
rownames(mds) = ss(rownames(mds), "-", 2)
pd = cbind(pd,mds[pd$IndividualID,])


#######
### all SNPs in locus "11", "129718630", "131718630
system(paste("plink --bfile", newbfile, 
	"--keep samples_to_extract.txt --snp rs10791097:130718630:T:G --window 1000",
	"--recode A-transpose --out all_snps_at_SNX19_GTEX --memory 225000"))
system(paste("plink --bfile", newbfile, 
	"--keep samples_to_extract.txt --snp rs10791097:130718630:T:G --window 1000",
	"--make-bed --out all_snps_at_SNX19_GTEX --memory 225000"))

## read in genotypes
genotypes2  = read.delim("all_snps_at_SNX19_DLPFC.traw",as.is=TRUE)
snpMapAll = genotypes2[,1:6]
snpMapAll$chrpos = paste0("chr", snpMapAll$CHR, ":", snpMapAll$POS)

genotypes2 = genotypes2[,-(1:6)]
rownames(genotypes2) = snpMapAll$SNP
colnames(genotypes2) = ss(colnames(genotypes2), "\\.",2)
snpAll = genotypes2[,pd$IndividualID]

save(pd, snpMat, snpInfo,snpAll, snpMapAll, compress=TRUE,
	file="GTEX_phenotype_annotated_SNX19_plusGenotype.rda")

	