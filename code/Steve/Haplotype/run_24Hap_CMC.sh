#$ -l mem_free=10G
#$ -l h_vmem=20G
#$ -l h_stack=256M
#$ -l h_fsize=1G
#$ -m e
#$ -M stephen.semick@gmail.com 

Rscript -e "rmarkdown::render('SNX19_HapStats_Analysis_24Hap_CMC.Rmd')"