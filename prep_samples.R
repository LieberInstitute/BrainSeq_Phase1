#####
require(GenomicRanges)

source("eqtl_functions.R")

# phenotype data
pd = read.csv("LIBD_RNA-seq_ALL_jhs_04_09_2014_toALL_750.csv", as.is=TRUE)[,1:33]
pd$RNum = paste0("R",pd$RNum)
pd$BRNum = paste0("Br",pd$BRNum)
colnames(pd)[2]= "BrNum"

# add bam info
pd$bamFile = paste0("/dcs01/lieber/ajaffe/Brain/DLPFC_PolyA/BAM/DLPFC_PolyA_",
	pd$RNum, "_accepted_hits.bam")
all(file.exists(pd$bamFile)) #  TRUE

## filter to 495 for the paper
keepIndex = which(pd$Dx %in% c("Control","Schizo") &
	pd$Race %in% c("AA","CAUC"))
pd = pd[keepIndex,]

############################
#### genotype data #########
fam = read.table("/dcs01/ajaffe/Imputation/Merged/imputed_plinkFiles/chr1.imputed.fam")
colnames(fam) = c("BrNum", "Platform", "MID","PID","Sex","Pheno")
pd = pd[which(pd$BrNum %in% fam$BrNum),] # keep samples w/ genotypes
fam = fam[!duplicated(fam$BrNum),] # remove 650 if they have 1M
famOut = fam[which(fam$BrNum %in% pd$BrNum),]
write.table(famOut[,1:2], "samples_to_extract.txt",
	col.names=FALSE, row.names=FALSE, quote=FALSE)

## merge by genotype,
# MAF > 5%, missingness < 10%, HWE p-value > 1e-6
bfile = "/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_imputed"
plinkPath = "/dcs01/ajaffe/Brain/DLPFC_PolyA/szControl/plink/"
system(paste0("plink --bfile ", bfile, " --keep samples_to_extract.txt --make-bed --geno 0.1 --maf 0.05  --hwe 0.000001 --out ", plinkPath, "LIBD_Brain_DLPFC_szControls_imputed --memory 225000"))

system(paste0("plink --bfile ", plinkPath, "LIBD_Brain_DLPFC_szControls_imputed", " --recode A-transpose --out ", plinkPath, "LIBD_Brain_DLPFC_szControls_imputed"))

## do MDS	
# ## ld prune from obs data
system(paste0("plink --bfile ", plinkPath, 
	"LIBD_Brain_DLPFC_szControls_imputed",
	" --indep 100 10 1.25 --geno 0.1 --out ", 
	plinkPath, "LIBD_Brain_DLPFC_szControls_imputed_indep"))

### MDS
system(paste0("plink --bfile ", plinkPath, 
	"LIBD_Brain_DLPFC_szControls_imputed --cluster --mds-plot 10 --extract ",
	plinkPath, "LIBD_Brain_DLPFC_szControls_imputed_indep.prune.in --out LIBD_szControl_n495"))

## read in genotypes
genotypes  = read.delim(paste0(plinkPath, 
	"LIBD_Brain_DLPFC_szControls_imputed.traw"), as.is=TRUE)
snpMap = genotypes[,1:6]
snpMap$CHR = gsub("23", "X", snpMap$CHR)

genotypes = genotypes[,-(1:6)]
rownames(genotypes) = snpMap$SNP
colnames(genotypes) = ss(colnames(genotypes), "_")
snp = genotypes[,pd$BrNum]

### add rs numbers
rs = read.table("/dcs01/ajaffe/Annotation/dbsnp142_common.txt",
	header = TRUE, comment.char="!", as.is = TRUE, sep = "\t")
rs$chromStart[rs$class=="single"] = rs$chromEnd[rs$class=="single"]

## overlap
nc = nchar(as.character(snpMap$ALT))
snpMapGR = GRanges(paste0("chr", snpMap$CHR), 
		IRanges(snpMap$POS, width=nc))
rsGR = GRanges(rs$chrom, IRanges(rs$chromStart, rs$chromEnd))	
oo = findOverlaps(snpMapGR, rsGR, type="equal")

## add rs
snpMap$name = NA
snpMap$name[queryHits(oo)] = rs$name[subjectHits(oo)]

## manually do others
ind = which(is.na(snpMap$name))
ind = ind[grep("^rs", snpMap$SNP[ind])]
snpMap$name[ind] = ss(as.character(snpMap$SNP[ind]), ":")

### add observed vs imputed
load("/users/ajaffe/Lieber/Projects/RNAseq/Consortium/PhaseI/observed_SNPs_byPlatform.rda")
obsList = lapply(infoList, function(x) paste0(x$chr, ":", x$position))
chrpos = paste0("chr",snpMap$CHR, ":", snpMap$POS)
byPlat = sapply(obsList, function(x) chrpos %in% x)
snpMap$numImp = 3-rowSums(byPlat)
snpMap$X.C.M = NULL

#####
mds = read.table("LIBD_szControl_n495.mds", 
	header=TRUE,as.is=TRUE,row.names=1)[,-(1:2)]
colnames(mds) = paste0("snpPC",1:ncol(mds))
pd = cbind(pd,mds[pd$BrNum,1:10])

## add total mapped
libSize = getTotalMapped(pd$bamFile,mc.cores=12)
pd$totalMapped = libSize$totalMapped
pd$mitoMapped = libSize$mitoMapped

## clean up variable names
pd$mitoRate = pd$mitoMapped / (pd$mitoMapped + pd$totalMapped)
names(pd)[c(5,17,18)] = c("Age", "totalAllMapped", "mappingRate")

##############
# bring in other demographic data
library(readxl)
pheno = read_excel("PhaseI_BrainSeq_Demos.xlsx",sheet=1)
colnames(pheno)[1] = "BrNum"
pheno$BrNum = paste0("Br", pheno$BrNum)
pheno = pheno[match(pd$BrNum, pheno$BrNum),]
colnames(pheno) = gsub(" ", "_", colnames(pheno))


## cause of death
pheno$Manner_Of_Death[pheno$Manner_Of_Death %in% c("No autopsy performed",
			"not filled in", "Pending", "Undetermined")] = NA
pd$Manner_Of_Death = pheno$Manner_Of_Death

# drugs
pheno$Antipsychotics[pheno$Antipsychotics %in% "Not Tested"] = NA
pd$Antipsychotics = pheno$Antipsychotics
pheno$"Antidepressants/_SSRIs"[pheno$"Antidepressants/_SSRIs" %in% "Not Tested"] = NA
pd$Antidepressants = pheno$"Antidepressants/_SSRIs"
pheno$Cotinine[pheno$Cotinine %in% c("not filled in", "Unknown")] = NA
pd$Cotinine = pheno$Cotinine
pheno$Nicotine[pheno$Nicotine %in% c("not filled in", "Unknown")] = NA
pd$Nicotine = pheno$Nicotine
pd$SmokingEither = ifelse(pd$Nicotine == "Positive" | 
	pd$Cotinine == "Positive", 	"Yes" ,"No")

# age of onset
pd$AgeOnsetSchizo = pheno$Age_Onset_Schizo

## PMI
pd$PMI = pheno$PMI

## save everything
save(pd, file="phenotype_annotated_szControlEqtl_DLPFC.rda")
save(snp, snpMap, file="/dcs01/ajaffe/Brain/DLPFC_PolyA/szControl/data/snpData_LIBD_szControls_cleaned.rda")

# write sample name
cat(pd$RNum, file="sample_rnums.txt",sep="\n")
