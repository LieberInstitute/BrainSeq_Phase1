##
library(ballgown)
library(rtracklayer)

ss =  ballgown:::ss
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/rdas/phenotype_annotated_szControlEqtl_DLPFC.rda")

## ballgown output from stringtie
gownDirs = paste0("/dcl01/lieber/ajaffe/Roche/CACNA1C/GTEX/StringTie/",
	pd$sra_accession)
table(file.exists(paste0(gownDirs, "/i2t.ctab"))) # TRUE

## read in
gownObj = ballgown(samples = gownDirs, 
	pData=pd, bamfiles=pd$bamFile)
save(gownObj, file=paste0(path, "ballgownObj_n495_full.rda"),
	compress=TRUE)
	
## filter
gownObjFilter = exprfilter(gownObj, cutoff=0.025)
save(gownObjFilter, file=paste0(path, "ballgownObj_n495_filter025.rda"),
	compress=TRUE)
	
#######
## extract and annotate
tFpkm = texpr(gownObjFilter, meas = "FPKM")
colnames(tFpkm) = ss(colnames(tFpkm), "\\.",2)
identical(pd$RNum, colnames(tFpkm))

## overlap to ensembl
tn = transcriptNames(gownObjFilter)
gtf = import(paste0(path, "/merged/merged.gtf"))
tInfo = mcols(gtf)[,c("transcript_id", "gene_id", 
	"gene_name", "nearest_ref", "class_code")]
tInfo = tInfo[!duplicated(tInfo$transcript_id),]
rownames(tInfo) = tInfo$transcript_id

round(100*prop.table(table(tInfo$class_code)),2)
    # =     i     j     o     p     r     s     u     x
# 32.04  0.02 56.23  1.15  0.00  0.02  0.01  9.09  1.43

## annotation
tMap = unlist(range(structure(gownObjFilter)$trans))
names(tMap) = rownames(tFpkm) = tn
mcols(tMap) = tInfo[match(names(tMap), rownames(tInfo)),]

round(100*prop.table(table(tMap$class_code)),2)
    # =     i     j     o     p     r     s     u     x
# 53.52  0.05 41.76  0.46  0.01  0.05  0.01  2.97  1.17

## add ensembl info
library(GenomicFeatures)
txdb = loadDb("../ensembl_v75_txdb.sqlite")
tx = transcriptsBy(txdb)
tx = unlist(tx)
tMap$EnsemblGeneID = names(tx)[match(tMap$nearest_ref, tx$tx_name)]

# save(tFpkm, tMap, file= "transcript_data_filtered_n495.rda")
save(tFpkm, tMap, file=paste0(path,"transcript_data_filtered_n495.rda"))
