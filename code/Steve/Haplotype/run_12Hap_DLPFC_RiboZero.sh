#$ -l mem_free=3G
#$ -l h_vmem=6G
#$ -l h_stack=256M
#$ -l h_fsize=1G
#$ -m e
#$ -M stephen.semick@gmail.com
#$ -o logs/$JOB_NAME.o
#$ -e logs/$JOB_NAME.e

Rscript -e "rmarkdown::render('SNX19_HapStats_Analysis_12Hap_DLPFC_RiboZero.Rmd')"