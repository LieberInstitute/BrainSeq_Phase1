#!/bin/sh

# qsub -t 294,299,315,482,486 -l mf=12G,h_vmem=16G,h_stack=256M rerun_stringtie_ballgown.sh
FILELIST=/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/samples_rnum.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
BAM=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/BAM/DLPFC_PolyA_${ID}_accepted_hits.bam
OUT=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/stringtie/$ID

# make outfile
mkdir -p $OUT

####### run stringtie ######
G=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/stringtie/merged/merged.gtf
stringtie $BAM -B -e -o $OUT/${ID}_merged_szControl -G $G