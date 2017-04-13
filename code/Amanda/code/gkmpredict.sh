# qsub -tc 245 -l mf=5G,h_vmem=5G,h_stack=256M -m e -M amanda.joy.price@gmail.com -t 1 gkmpredict.sh

cd /dcl01/lieber/ajaffe/Amanda/SNX19/code/Amanda/
/users/aprice26/biotools/lsgkm-master/bin/gkmpredict -T 16 supporting_files/ref_extended_Hap.fa supporting_files/neuronsTogether_replicated.300bp.pval5_repeat.model.txt data/refSNPs_11mer_gkmpredict.txt
/users/aprice26/biotools/lsgkm-master/bin/gkmpredict -T 16 supporting_files/alt_extended_Hap.fa supporting_files/neuronsTogether_replicated.300bp.pval5_repeat.model.txt data/altSNPs_11mer_gkmpredict.txt
