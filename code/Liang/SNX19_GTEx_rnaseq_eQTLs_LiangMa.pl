#!/usr/bin/perl -w

use strict;


my $f_data = "rawAndRpkmCounts_plusGenotype_SNX19_GTEX_n8206.rda";
my $d_bat  = "";
my $d_res  = "";

my $f_list = "tissue_names.txt";



my $cmd;

open LIST, $f_list or die "Cannot open $f_list";
while(<LIST>) {

		chomp;
		chop;
		
		my $tissue = $_;
		
		my $tissue_name = $tissue;
		
		$tissue_name =~ s/\s+//g;
		$tissue_name =~ s/\(/_/g;
		$tissue_name =~ s/\)//g;
		
    my $f_bat = "$d_bat/run_tissue_$tissue_name.R";
    my $f_out = "$d_res/run_tissue_$tissue_name.RData";
    my $f_job = "$d_bat/run_tissue_$tissue_name.sh";

    open(BAT, ">$f_bat") or die("cannot open $f_bat");
    
    print BAT "rm(list=ls())\n";

    print BAT "library(MatrixEQTL)\n";
    print BAT "library(GenomicRanges)\n";
    print BAT "library(sva)\n";
    print BAT "library(XVector)\n";
    print BAT "library(Biostrings)\n";

    print BAT "load(\"$f_data\")\n";

    print BAT "bgGene = matrix(rep(pd\$totalMapped), nc = nrow(pd), nr = nrow(geneCounts),   byrow=TRUE)\n";
    print BAT "widGGene = matrix(rep(geneMap\$Length), nr = nrow(geneCounts), nc = nrow(pd), byrow=FALSE)\n";
    print BAT "geneRpkm = geneCounts/(widGGene/1000)/(bgGene/1e6)\n";

    print BAT "bgJunc = matrix(rep(pd\$totalMapped/80e6), nc = nrow(pd), nr = nrow(jCounts), byrow=TRUE)\n";
    print BAT "jRp80M = jCounts/bgJunc\n";

    print BAT "tissues = unique(pd\$SMTSD)\n";
    print BAT "tissue = tissues[grep(\"$tissue\", tissues)]; \n";
 
    print BAT "aIndex = which((pd\$AGE>13) & (pd\$SMTSD==tissue))\n";
    print BAT "pd2 = pd[aIndex,]\n";
    print BAT "snp2 = as.matrix(snpAll[,aIndex])\n";
    print BAT "geneRpkm2 = as.matrix(log2(geneRpkm[,aIndex]+1))\n";
    print BAT "exonRpkm2 = as.matrix(log2(exonRpkm[,aIndex]+1))\n";
    print BAT "jRp80M2 = as.matrix(log2(jRp80M[,aIndex]+1))\n";
    
    print BAT "exprs2 = rbind(geneRpkm2, exonRpkm2, jRp80M2);\n";
    
    print BAT "geneCounts2 = geneCounts[,aIndex]\n";
    print BAT "exonCounts2 = exonCounts[,aIndex]\n";
    print BAT "jCounts2 = jCounts[,aIndex]\n";
    
    print BAT "geneCounts_avg = rowMeans(geneCounts2)\n";
    print BAT "exonCounts_avg = rowMeans(exonCounts2)\n";
    print BAT "jCounts_avg    = rowMeans(jCounts2)\n";
    print BAT "avg_counts = c(geneCounts_avg, exonCounts_avg, jCounts_avg);\n";
    print BAT "avg_counts = data.frame(avg_counts)\n";
    print BAT "avg_counts\$Feature = rownames(avg_counts)\n";
    
    print BAT "gm = cbind(rownames(geneMap), geneMap\$Symbol)\n";
    print BAT "colnames(gm) = c(\"GeneID\", \"Symbol\")\n";
    print BAT "em = exonMap[,c(\"Geneid\", \"Symbol\")]\n";
    print BAT "colnames(em)[1] = \"GeneID\"\n";
    print BAT "jm = cbind(jMap\$newGeneID, jMap\$newGeneSymbol)\n";
    print BAT "colnames(jm) = c(\"GeneID\", \"Symbol\")\n";
    
    print BAT "jointMap = rbind(gm,em,jm)\n";

    print BAT "jointMap\$GeneID = as.character(jointMap\$GeneID)\n";
    print BAT "jointMap\$Symbol = as.character(jointMap\$Symbol)\n";
    print BAT "jointMap\$Feature = rownames(exprs2)\n";
    print BAT "dim(jointMap) # 8044    3\n";
    
    print BAT "jointMap\$featureType = rep(c(\"Gene\", \"Exon\", \"Jxn\"), times = c(nrow(geneRpkm), nrow(exonRpkm),nrow(jRpkm))) \n";  
    print BAT "dim(jointMap) \n";

    print BAT "jointMap2 = merge(jointMap, avg_counts, by=\"Feature\", sort=F) \n";

    print BAT "pca = prcomp(t(exprs2)) \n";
    print BAT "exprsPC = pca\$x[,1:15] \n";

    print BAT "mod = model.matrix(~pd2\$snpPC1 + pd2\$snpPC2 + pd2\$snpPC3 + pd2\$snpPC4 + pd2\$snpPC5 + exprsPC)  \n";
    print BAT "colnames(mod)[2:ncol(mod)] = c(paste0(\"snpPC\",1:5), paste0(\"exprsPC\",1:15)) \n";

    print BAT "covs = SlicedData\$new(t(mod[,-1])) ##delete the first colume: intercect. Left snpPC1-5 and exonPCs1-15. \n";

    print BAT "snpspos = snpMapAll[,c(\"SNP\",\"CHR\",\"POS\")] #the snps name in output defined here!! \n";
    print BAT "snpspos\$CHR = paste0(\"chr\",snpspos\$CHR) # past 'chr' to the 'CHR' colume. \n";
    print BAT "colnames(snpspos) = c(\"name\",\"chr\",\"pos\") \n";

    print BAT "posGene = geneMap[,1:3]\n";
    print BAT "posGene\$name = rownames(geneMap)\n";
    print BAT "posGene = posGene[,c(4,1:3)]\n";

    print BAT "posExon = exonMap[,2:4]\n";
    print BAT "posExon\$name = rownames(exonMap)\n";
    print BAT "posExon = posExon[,c(4,1:3)]\n";

    print BAT "posJxn = as.data.frame(jMap)[,1:3]\n";
    print BAT "posJxn\$name = names(jMap)\n";
    print BAT "posJxn = posJxn[,c(4,1:3)]\n";
    print BAT "names(posJxn)[2:4] = c(\"Chr\", \"Start\",\"End\")\n";

    print BAT "snpSlice = SlicedData\$new(snp2) # SNP 2831\n";
    print BAT "snpSlice\$ResliceCombined(sliceSize = 5000)\n";

    print BAT "exprsSlice = SlicedData\$new(exprs2) \n";
    print BAT "pos = rbind(posGene, posExon, posJxn)\n";

    print BAT "meJoint = Matrix_eQTL_main(snps=snpSlice, gene = exprsSlice,\n";
    print BAT "    cvrt = covs, output_file_name.cis =  \".ctxt\" ,\n";
    print BAT "    pvOutputThreshold.cis = 1,  pvOutputThreshold=0,\n";
    print BAT "    snpspos = snpspos, genepos = pos, \n";
    print BAT "    useModel = modelLINEAR,  cisDist=5e6,\n";
    print BAT "    pvalue.hist = 100, min.pv.by.genesnp = T)\n";

    print BAT "all = meJoint\$cis\$eqtls\n";

    print BAT "coefMat = merge(jointMap2, all, by.x=\"Feature\", by.y=\"gene\")\n";

    print BAT "coefMat\$Tissue = tissue\n";
    print BAT "coefMat\$SampleSize = length(aIndex)\n";
    print BAT "save(coefMat, file=\"$f_out\")\n";

    close(BAT);

    open JOB, ">$f_job" or die "Cannot open $f_job!!";
    print JOB "R CMD BATCH $f_bat\n";
    close(JOB);

    $cmd = "qsub -cwd -l mem_free=80G,h_vmem=100G $f_job";
    system($cmd);


    # print "R CMD BATCH $f_bat\n";
    # system("R CMD BATCH $f_bat");


}

close LIST;


