##
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
rm(geneRpkm, geneMap, exonMap,exonRpkm)

## only controls 
keepIndex = which(pd$Dx == "Control")
pd = pd[keepIndex,]
jRpkm = jRpkm[,keepIndex]

## drop all 0s 
exprsIndex = which(rowSums(jRpkm) > 0)
jRpkm = jRpkm[exprsIndex,]
jMap = jMap[exprsIndex]

# groups
pd$metricGroups = "Adult"
pd$metricGroups[pd$Age < 0] = "Fetal"
pd$metricGroups[pd$Age > 0 & pd$Age < 13] = "Child"
pd$metricGroups = factor(pd$metricGroups, 
	levels = c("Fetal","Child","Adult"))

###############
## junctions 
###############
library(rtracklayer)

## means 
gIndexes = split(1:nrow(pd), pd$metricGroups)
jMeans = sapply(gIndexes, function(ii) rowMeans(jRpkm[,ii]))

## covert means to scores
jMeansRound = round(jMeans)
jMeansRound[jMeansRound >= 1000] = 999

## write out BEDs 	
tmp = as.data.frame(jMap) 
tmp$score = tmp$bedColorsRgb = NA
tmp$bedSign = ifelse(jMap$code == "InEns", "+", "-")
tmp$name = rownames(tmp)

bed = tmp[,c("seqnames", "start", "end", 
	"name", "score", "bedSign", "start", "end",
	"bedColorsRgb")]
colnames(bed) = c("chr","start","end","name","score", "strand",
		"thickStart", "thickEnd", "itemRgb")

header = paste0("track name=LIBD_jxn_", colnames(jMeans),
	" description='Junctions from LIBD RNA-seq - ", colnames(jMeans),
	"' visibility=2 itemRgb='On'")
fn = paste0("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/hub/hg19/",
		colnames(jMeans), "_junctions.bed")		
		
## write out each file
for(i in 1:ncol(jMeansRound)) {
	cat(".")
	out = bed
	out$score = jMeansRound[,i]
	
	cols = colorRampPalette(c("white","red"))(1000)[out$score+1]
	cols[bed$strand == "+"] = colorRampPalette(c("white","black"))(1000)[out$score+1][bed$strand == "+"]
	cols = t(col2rgb(cols))
	cols = paste(cols[,1], cols[,2], cols[,3], sep=",")
	out$itemRgb = cols
	
	# no header for big bed	
	# write.table(header[i], file=fn[i],	row.names=FALSE, 
		# col.names=FALSE, quote=FALSE,append=TRUE)
	write.table(out, file=fn[i], row.names=FALSE, 
		col.names=FALSE, quote=FALSE,append=FALSE)	
}

### convert to bigbed
# system("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/fetchChromSizes hg19 > hg19.chrom.sizes")
o = gsub(".bed",".bigBed", fn,fixed=TRUE)
thecall = paste("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/bedToBigBed -type=bed9",
	fn, "hg19.chrom.sizes", o) 
sapply(thecall, system)

################
## coverage ####
library(derfinder)

## normalization factor
pd$mappedPer80M = pd$totalMapped/80e6

## load full coverage
load("/dcl01/lieber/ajaffe/derRuns/controlEQTL/control/CoverageInfo/fullCov.Rdata")
fullCov = lapply(fullCov, function(y) y[,pd$RNum]) # order

## take mean
meanCov = lapply(fullCov, function(y) {
	cat(".")
	yM = DataFrame(mapply(function(x,d) x/d, y, pd$mappedPer80M))
	means = sapply(gIndexes, function(ii) Reduce("+", yM[ii]) / length(ii))
})

# export
for(i in seq(along=gIndexes)) {
	cat(names(gIndexes)[i])
	gcov = lapply(meanCov, "[[", i)
	gcov = lapply(gcov, function(x) {
		runValue(x) = signif(runValue(x),3)
		x})
	out = RleList(gcov,compress=FALSE)
	export(out, con=paste0("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/hub/hg19/",
		colnames(jMeans)[i], "_coverage.bw"))
	cat("\n")
}

system("nano /dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/hub/hg19/genomes.txt")
# genome hg19
# trackDb hg19/trackDb.txt

system("nano /dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/hub/hg19/hub.txt")
# hub DLPFC_RNAseq 
# shortLabel LIBD Human DLPFC RNA-seq Expression
# longLabel Junctions and coverage from RNAseq data across human brain development
# genomesFile genomes.txt
# email andrew.jaffe@libd.org
# descriptionUrl http://www.libd.org
