###
# qsub -tc 175 -l mf=4G,h_vmem=6G,h_stack=256M -t 1-495 run_degradation_coverage.sh
FILELIST=/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/sample_rnums.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
BED=/users/ajaffe/Lieber/Projects/RNAseq/Degradation/regionApproach/bed/polyA_degradation_regions_v2.bed
BAM=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/BAM/DLPFC_PolyA_${ID}_accepted_hits.bam
BW=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/Coverage/DLPFC_PolyA_${ID}.bw
BG=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/Coverage/DLPFC_PolyA_${ID}.bedGraph
CHRSIZE=/users/ajaffe/Lieber/Projects/RNAseq/hg19.chrom.sizes

if [ ! -e $BW ] ; then
	bedtools genomecov -ibam ${BAM} -bga -split > $BG
	bedGraphToBigWig $BG $CHRSIZE $BW
	rm $BG
fi 

OUT=/users/ajaffe/Lieber/Projects/RNAseq/Consortium/PhaseI/degradation/DLPFC_PolyA_R${ID}_degradeStats.txt
bwtool summary $BED $BW $OUT -header -fill=0 -with-sum
