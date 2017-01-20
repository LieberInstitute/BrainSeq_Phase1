###
library(GenomicRanges)
library(limma)
source("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl_functions.R")

#################################
#### LOAD DATA ##################
#################################

# load SNPs for annotation info
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/snpData_LIBD_szControls_cleaned.rda")
snpMap$CHR[snpMap$CHR == "23"] = "X"

### add number minor alleles
snpMap$numMinorHom = rowSums(snp == 2,na.rm=TRUE)
snpMap$numMinorHom_13plus = rowSums(snp[,pd$Age > 13] == 2,na.rm=TRUE)

# load expression annotation
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
exonMap$coord = with(exonMap, paste0(Chr, ":", Start, "-", End, "(", Strand, ")"))

# ERs + Tx
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/transcript/transcript_data_filtered_n495.rda")

oo = findOverlaps(regions, makeGRangesFromDataFrame(geneMap))
regions$EnsemblID = NA
regions$EnsemblID[queryHits(oo) ] = rownames(geneMap)[subjectHits(oo)]

################################
## bring back in eqtls #########
################################

###### GENE
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/annotated_gene_eqtl_szControl_cisOnly.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/gene_eqtl_szControl_cisOnly.rda")
geneEqtl = meGeneCis$cis$eqtls
geneEqtl = geneEqtl[geneEqtl$FDR < 0.01,]
geneEqtl$gene = as.character(geneEqtl$gene)
geneEqtl$snps = as.character(geneEqtl$snps)
geneEqtl$bonf = p.adjust(geneEqtl$pvalue, "bonf", 
	n = (sigGene$bonf[1]/sigGene$pvalue[1]))

###### EXON
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/annotated_exon_eqtl_szControl_cisOnly.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/exon_eqtl_szControl_cisOnly.rda")
exonEqtl = meExonCis$cis$eqtls
exonEqtl = exonEqtl[exonEqtl$FDR < 0.01,]
exonEqtl$gene = as.character(exonEqtl$gene)
exonEqtl$snps = as.character(exonEqtl$snps)
exonEqtl$bonf = p.adjust(exonEqtl$pvalue, "bonf", 
	n = (sigExon$bonf[1]/sigExon$pvalue[1]))

######### TX
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/annotated_transcript_eqtl_szControl_cisOnly.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/transcript_eqtl_szControl_cisOnly.rda")
transcriptEqtl = meTransCis$cis$eqtls
transcriptEqtl = transcriptEqtl[transcriptEqtl$FDR < 0.01,]
transcriptEqtl$gene = as.character(transcriptEqtl$gene)
transcriptEqtl$snps = as.character(transcriptEqtl$snps)
transcriptEqtl$bonf = p.adjust(transcriptEqtl$pvalue, "bonf", 
	n = (sigTrans$bonf[1]/sigTrans$pvalue[1]))
	
###### JUNCTION
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/annotated_junction_eqtl_szControl_cisOnly.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/junction_eqtl_szControl_cisOnly.rda")
junctionEqtl = meJxnCis$cis$eqtls
junctionEqtl = junctionEqtl[junctionEqtl$FDR < 0.01,]
junctionEqtl$gene = as.character(junctionEqtl$gene)
junctionEqtl$snps = as.character(junctionEqtl$snps)
junctionEqtl$bonf = p.adjust(junctionEqtl$pvalue, "bonf", 
	n = (sigJxn$bonf[1]/sigJxn$pvalue[1]))

###### DER
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/annotated_der_eqtl_szControl_cisOnly.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/der_eqtl_szControl_cisOnly.rda")
derEqtl = meDerCis$cis$eqtls
derEqtl = derEqtl[derEqtl$FDR < 0.01,]
derEqtl$gene = as.character(derEqtl$gene)
derEqtl$snps = as.character(derEqtl$snps)
derEqtl$bonf = p.adjust(derEqtl$pvalue, "bonf", 
	n = (sigDer$bonf[1]/sigDer$pvalue[1]))

### add number of homozygous minor
geneEqtl$numMinor13 = snpMap$numMinorHom_13plus[
	match(geneEqtl$snps, snpMap$SNP)]
exonEqtl$numMinor13 = snpMap$numMinorHom_13plus[
	match(exonEqtl$snps, snpMap$SNP)]
transcriptEqtl$numMinor13 = snpMap$numMinorHom_13plus[
	match(transcriptEqtl$snps, snpMap$SNP)]
junctionEqtl$numMinor13 = snpMap$numMinorHom_13plus[
	match(junctionEqtl$snps, snpMap$SNP)]
derEqtl$numMinor13 = snpMap$numMinorHom_13plus[
	match(derEqtl$snps, snpMap$SNP)]

###########################
### add gene and transcript info
###############################

geneEqtl$Symbol = geneMap$Symbol[match(geneEqtl$gene, rownames(geneMap))]
geneEqtl$EnsemblID = as.character(geneEqtl$gene)
geneEqtl$Class = "InEns"

exonEqtl$Symbol = exonMap$Symbol[match(exonEqtl$gene, rownames(exonMap))]
exonEqtl$EnsemblID = exonMap$Geneid[match(exonEqtl$gene, rownames(exonMap))]
exonEqtl$Class = "InEns"

transcriptEqtl$Symbol = tMap$gene_name[match(transcriptEqtl$gene, names(tMap))]
transcriptEqtl$EnsemblID = tMap$EnsemblGeneID[match(transcriptEqtl$gene, names(tMap))]
code = tMap$class_code[match(transcriptEqtl$gene, names(tMap))]
transcriptEqtl$Class = "Novel"
transcriptEqtl$Class[code == "="] = "InEns"
transcriptEqtl$Class[code == "j"] = "ExonSkip"
transcriptEqtl$Class[code %in% c("c", "e", "o", "p", "x") ] = "AltStartEnd"

junctionEqtl$Symbol = jMap$newGeneSymbol[match(junctionEqtl$gene, names(jMap))]
junctionEqtl$EnsemblID = jMap$newGeneID[match(junctionEqtl$gene, names(jMap))]
junctionEqtl$Class = jMap$code[match(junctionEqtl$gene, names(jMap))]

derEqtl$Symbol = regions$nearestSymbol[match(derEqtl$gene, names(regions))]
derEqtl$EnsemblID = regions$EnsemblID[match(derEqtl$gene, names(regions))]
code2 = regions$annoClass[match(derEqtl$gene, names(regions))]
derEqtl$Class = "Novel"
derEqtl$Class[code2 == "strictExonic"] = "InEns"
derEqtl$Class[code2 %in% c("exonIntron", "extendUTR")] = "AltStartEnd"

###################
### try merging ###
allEqtls = rbind(geneEqtl, exonEqtl, transcriptEqtl, junctionEqtl, derEqtl)
allEqtls$snpRsNum = snpMap$name[match(allEqtls$snps, snpMap$SNP)]
allEqtls$Type = rep(c("Gene","Exon", "Transcript", "Junction", "ER"),
	c(nrow(geneEqtl), nrow(exonEqtl), nrow(transcriptEqtl), 
		nrow(junctionEqtl), nrow(derEqtl)))
allEqtls$statistic = NULL

## make smaller 
allEqtls = DataFrame(allEqtls)
allEqtls$Class = Rle(allEqtls$Class)
allEqtls$Type = Rle(allEqtls$Type)
allEqtls$Type = Rle(allEqtls$Type)
allEqtls = allEqtls[order(allEqtls$pvalue),]
names(allEqtls)[1:2] = c("SNP", "Feature")

#############################
### get transcript info #####
xx=load("/users/ajaffe/Lieber/Projects/RNAseq/ensembl_V75_feature_to_Tx.rda")

## add annotation from tx
theTxs = CharacterList(as.list(tMap$nearest_ref))
names(theTxs) = names(tMap)
allTx = c(allTx, theTxs)

mmTx = match(allEqtls$Feature, names(allTx))
tx = CharacterList(vector("list", nrow(allEqtls)))
tx[!is.na(mmTx)] = allTx[mmTx[!is.na(mmTx)]]
allEqtls$NumTxEqtl = elementLengths(tx)
allEqtls$WhichTx = tx

mmTxGene = match(allEqtls$EnsemblID, names(allTx))
txGene = CharacterList(vector("list", nrow(allEqtls)))
txGene[!is.na(mmTxGene)] = allTx[mmTxGene[!is.na(mmTxGene)]]
allEqtls$NumTxGene = elementLengths(txGene)

allEqtls$NumTxEqtl = Rle(allEqtls$NumTxEqtl)
allEqtls$NumTxGene = Rle(allEqtls$NumTxGene)

###################################
########## bonf filter ############

bonfEqtls = allEqtls[which(allEqtls$bonf < 0.05),]
bonfEqtls$Type = as.character(bonfEqtls$Type)

###########################
#### map to unique SNPs ###
uSnps = unique(bonfEqtls$SNP)

uniqueSnps = as.data.frame(matrix(FALSE, nr = length(uSnps), nc = 5,
	dimnames = list(uSnps, c("Gene","Exon", "Transcript", "Junction", "ER"))))
uniqueSnps$Gene = rownames(uniqueSnps) %in% bonfEqtls$SNP[bonfEqtls$Type=="Gene"]
uniqueSnps$Exon = rownames(uniqueSnps) %in% bonfEqtls$SNP[bonfEqtls$Type=="Exon"]
uniqueSnps$Transcript = rownames(uniqueSnps) %in% bonfEqtls$SNP[bonfEqtls$Type=="Transcript"]
uniqueSnps$Junction = rownames(uniqueSnps) %in% bonfEqtls$SNP[bonfEqtls$Type=="Junction"]
uniqueSnps$ER = rownames(uniqueSnps) %in% bonfEqtls$SNP[bonfEqtls$Type=="ER"]

colMeans(uniqueSnps)
table(rowSums(uniqueSnps))

###########################
#### map to unique Tx ###
uTxs = unique(unlist(bonfEqtls$WhichTx))

uniqueTxs = as.data.frame(matrix(FALSE, nr = length(uTxs), nc = 5,
	dimnames = list(uTxs, c("Gene","Exon", "Transcript", "Junction", "ER"))))
uniqueTxs[unlist(bonfEqtls$WhichTx[bonfEqtls$Type == "Gene"]),"Gene"] = TRUE
uniqueTxs[unlist(bonfEqtls$WhichTx[bonfEqtls$Type == "Exon"]),"Exon"] = TRUE
tt1 = unlist(bonfEqtls$WhichTx[bonfEqtls$Type == "Transcript"])
tt1 = tt1[!is.na(tt1)]
uniqueTxs[tt1,"Transcript"] = TRUE
uniqueTxs[unlist(bonfEqtls$WhichTx[bonfEqtls$Type == "Junction"]),"Junction"] = TRUE
uniqueTxs[unlist(bonfEqtls$WhichTx[bonfEqtls$Type == "ER"]),"ER"] = TRUE

colMeans(uniqueTxs)
table(rowSums(uniqueTxs))

###########################
#### map to unique Gene ###
uGenes = unique(bonfEqtls$EnsemblID)
uniqueGenes = as.data.frame(matrix(FALSE, nr = length(uGenes), nc = 5,
	dimnames = list(uGenes, c("Gene","Exon", "Transcript", "Junction", "ER"))))
uniqueGenes$Gene = rownames(uniqueGenes) %in% bonfEqtls$EnsemblID[bonfEqtls$Type=="Gene"]
uniqueGenes$Exon = rownames(uniqueGenes) %in% bonfEqtls$EnsemblID[bonfEqtls$Type=="Exon"]
uniqueGenes$Transcript = rownames(uniqueGenes) %in% bonfEqtls$EnsemblID[bonfEqtls$Type=="Transcript"]
uniqueGenes$Junction = rownames(uniqueGenes) %in% bonfEqtls$EnsemblID[bonfEqtls$Type=="Junction"]
uniqueGenes$ER = rownames(uniqueGenes) %in% bonfEqtls$EnsemblID[bonfEqtls$Type=="ER"]
colMeans(uniqueGenes)
table(rowSums(uniqueGenes))

colSums(uniqueGenes[rowSums(uniqueGenes) == 1,])

##########################
#### snp - gene ###########
bonfEqtls$snpGene = paste0(bonfEqtls$SNP, ".", bonfEqtls$EnsemblID)

uSnpGenes = unique(bonfEqtls$snpGene)
uniqueSnpGenes = as.data.frame(matrix(FALSE, nr = length(uSnpGenes), nc = 5,
	dimnames = list(uSnpGenes, c("Gene","Exon", "Transcript", "Junction", "ER"))))

uniqueSnpGenes$Gene = rownames(uniqueSnpGenes) %in% bonfEqtls$snpGene[bonfEqtls$Type=="Gene"]
uniqueSnpGenes$Exon = rownames(uniqueSnpGenes) %in% bonfEqtls$snpGene[bonfEqtls$Type=="Exon"]
uniqueSnpGenes$Transcript = rownames(uniqueSnpGenes) %in% bonfEqtls$snpGene[bonfEqtls$Type=="Transcript"]
uniqueSnpGenes$Junction = rownames(uniqueSnpGenes) %in% bonfEqtls$snpGene[bonfEqtls$Type=="Junction"]
uniqueSnpGenes$ER = rownames(uniqueSnpGenes) %in% bonfEqtls$snpGene[bonfEqtls$Type=="ER"]

colMeans(uniqueSnpGenes)
table(rowSums(uniqueSnpGenes))

uList = split(uniqueSnpGenes, factor(ss(rownames(uniqueSnpGenes),
	"\\.",2),	levels = unique(bonfEqtls$EnsemblID)))
uToGene = t(sapply(uList, function(x) table(factor(rowSums(x),levels=1:5))))
uList3 = split(uniqueSnpGenes[,c(1,2,4)], factor(ss(rownames(uniqueSnpGenes),
	"\\.",2),	levels = unique(bonfEqtls$EnsemblID)))
uToGene3 = t(sapply(uList3, function(x) table(factor(rowSums(x),levels=0:3))))

## clean expression, filter for age
aIndex= which(pd$Age > 13)
pd2= pd[aIndex,]
snp2 = as.matrix(snp[,aIndex])
geneRpkm2 = as.matrix(log2(geneRpkm[,aIndex]+1))
exonRpkm2 = as.matrix(log2(exonRpkm[,aIndex]+1))
jRpkm2 = as.matrix(log2(jRpkm[,aIndex]+1))
regionMat2 = as.matrix(log2(regionMat[,aIndex]+1))
tFpkm2 = as.matrix(log2(tFpkm[,aIndex]+1))

## statistical model
mod = model.matrix(~pd2$snpPC1 + pd2$snpPC2 + 
	pd2$snpPC3 + pd2$snpPC4 + pd2$snpPC5)