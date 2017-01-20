##

# qsub -tc 120 -l mf=4G,h_vmem=6G,h_stack=256M -t 1-495 counts_junctions.sh
FILELIST=/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/sample_rnums.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
BAM=/dcs01/lieber/ajaffe/Brain/DLPFC_PolyA/BAM/DLPFC_PolyA_${ID}_accepted_hits.bam
OUTJXN=/dcs01/lieber/ajaffe/Brain/DLPFC_PolyA/Junctions/DLPFC_PolyA_${ID}_junctions_primaryOnly.count
SCRIPT=/dcl01/leek/data/introns_from_sam.py
PYPY=/home/student/anellor1/raildotbio/pypy-2.5-linux_x86_64-portable/bin/pypy
samtools view -bh -F 0x100 $BAM | $PYPY $SCRIPT > $OUTJXN