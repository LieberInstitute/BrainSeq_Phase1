##

# qsub -pe local 3 -t 1-495 -tc 50 -l mf=8G,h_vmem=12G,h_stack=256M run_stringtie.sh
# sed -e 's/^/chr/' /dcl01/lieber/ajaffe/Annotation/Homo_sapiens.GRCh37.75.gtf > /dcl01/lieber/ajaffe/Annotation/Homo_sapiens.GRCh37.75_chrPrefix.gtf

FILELIST=/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/samples_rnum.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST)
BAM=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/BAM/DLPFC_PolyA_${ID}_accepted_hits.bam

OUT=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/stringtie/$ID
mkdir -p $OUT
a=/dcl01/lieber/ajaffe/Annotation/Homo_sapiens.GRCh37.75_chrPrefix.gtf
stringtie -p 3 $BAM -o $OUT/${ID}_stringtie_ensemblV75.gtf -G $a 