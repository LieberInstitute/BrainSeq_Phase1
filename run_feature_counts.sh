###

# qsub -tc 120 -l mf=4G,h_vmem=6G,h_stack=256M -t 1-495 run_feature_counts.sh
FILELIST=/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/sample_rnums.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
BAM=/dcs01/lieber/ajaffe/Brain/DLPFC_PolyA/BAM/DLPFC_PolyA_${ID}_accepted_hits.bam

featureCounts -A /dcl01/lieber/ajaffe/Annotation/chrAliases_GRCh37_to_hg19.csv \
	-a /dcl01/lieber/ajaffe/Annotation/Homo_sapiens.GRCh37.75.gtf \
	-o /dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/Counts/Gene/DLPFC_PolyA_${ID}_Ensembl_v75_Genes.counts $BAM
featureCounts -O -f -A /dcl01/lieber/ajaffe/Annotation/chrAliases_GRCh37_to_hg19.csv \
	-a /dcl01/lieber/ajaffe/Annotation/Homo_sapiens.GRCh37.75.gtf \
	-o /dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/Counts/Exon/DLPFC_PolyA_${ID}_Ensembl_v75_Exons.counts $BAM
	