find /dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/stringtie -name R\*.gtf -print > all_gtfs.txt
a=/dcl01/lieber/ajaffe/Annotation/Homo_sapiens.GRCh37.75_chrPrefix.gtf
out=/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/stringtie/merged
chr=/dcl01/lieber/ajaffe/Annotation/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/
cuffmerge -p 4 -o $out -g $a -s $chr all_gtfs.txt
