##

# qsub -tc 150 -l mf=10G,h_vmem=15G,h_stack=256M -t 1-251 count_junctions_regtool.sh
FILELIST=/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/rerun.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
BAM=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/BAM/DLPFC_PolyA_${ID}_accepted_hits.bam
OUTJXN=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/Junctions/DLPFC_PolyA_${ID}_junctions_primaryOnly_regtool.bed
TMPBAM=$TMPDIR/${ID}.bam
samtools view -bh -F 0x100 $BAM > $TMPBAM
samtools index $TMPBAM
regtools junctions extract -i 9 -o $OUTJXN $TMPBAM
OUTCOUNT=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/Junctions/DLPFC_PolyA_${ID}_junctions_primaryOnly_regtool.count
/users/ajaffe/Lieber/Projects/RNAseq/bed_to_juncs_withCount < $OUTJXN > $OUTCOUNT