####
source("../../eqtl_functions.R")

# load output
load("rdas/annotated_gene_eqtl_control_13plus_cisOnly.rda")
load("rdas/annotated_exon_eqtl_control_13plus_cisOnly.rda")
load("rdas/annotated_junction_eqtl_control_13plus_cisOnly.rda")
load("rdas/annotated_transcript_eqtl_control_13plus_cisOnly.rda")
load("rdas/annotated_der_eqtl_control_13plus_cisOnly.rda")

# unique SNPs?
theSnps = c(sigGene$snps, sigExon$snps, sigDer$snps,
	sigTrans$snps, sigJxn$snps)

# make list
eqtlList = list(Gene = sigGene, Exon = sigExon, 
	Transcript =sigTrans, Junction = sigJxn, DER = sigDer)

# FDR
sapply(eqtlList, function(x) sum(x$FDR < 0.01))	
# bonf
sapply(eqtlList, function(x) sum(x$bonf < 0.05))	

## significant ones at bonf
sigList = lapply(eqtlList, 
	function(x) x[which(x$bonf < 0.05),])

# p-value cutoff for BONF < 5%
sapply(eqtlList, function(x) max(x$pvalue))
sapply(sigList, function(x) max(x$pvalue))

# number of unique genes
length(unique(sigList$Exon$exprs_geneid))
length(unique(sigList$Transcript$exprs_EnsemblGeneID[
	!is.na(sigList$Transcript$exprs_EnsemblGeneID)]))
length(unique(sigList$Junction$exprs_newGeneID))
length(unique(sigList$DER$exprs_symbol[abs(sigList$DER$exprs_distance) == 0]))

# mean maf
sapply(sigList, function(x) mean(x$inSampleMAF))

## effect size
t(sapply(sigList, function(x) 2^quantile(abs(x$beta))[c(3,2,4)]))


## mean expression
sapply(sigList, function(x) mean(x[,grep("mean", names(x))]))

# dist
t(sapply(sigList, function(x) quantile(x[,grep("istTo", names(x))],
	na.rm=TRUE)[c(3,2,4)]))

## SNPs in LD
t(sapply(sigList, function(x) quantile(x$numSnps)[c(3,2,4)]))

## SNP range
t(sapply(sigList, function(x) 
	round(quantile(x$snpLDLength)[c(3,2,4)]/1e3)))

########novel
## transcript
round(prop.table(table(sigList$Transcript$exprs_class_code))*100,2)

## junction
prop.table(table(sigList$Junction$exprs_code))

## der
mean(sigList$DER$exprs_exon > 0)
mean(sigList$DER$exprs_exon == 0)

########################
#### transcript specific
tList = split(sigList$Transcript$exprs_nearest_ref, 
	sigList$Transcript$exprs_EnsemblGeneID)
numTrans = table(sapply(tList,length))
numTrans[1]/ sum(numTrans)

## junction
prop.table(table(sigList$Junction$exprs_code))

table(sigList$Junction$exprs_numTx,useNA="ifany")
round(prop.table(table(sigList$Junction$exprs_numTx,useNA="ifany")),3)
sum(table(sigList$Junction$exprs_numTx)[1])/nrow(sigList$Junction)
sum(table(sigList$Junction$exprs_numTx)[1:2])/nrow(sigList$Junction)

#################### 
## all genes ####
allGenesSym = unique(c(sigList$Gene$exprs_symbol, sigList$Exon$exprs_symbol,
	sigList$Junction$exprs_newGeneSym, sigList$Transcript$exprs_ref_gene_id,
	sigList$DER$exprs_symbol[sigDer$exprs_distance == 0]))
allGenesSym = allGenesSym[!is.na(allGenesSym) & !grepl("-", allGenesSym)]
length(allGenesSym)

allGenesID = unique(c(sigList$Gene$gene, sigList$Exon$exprs_geneid,
	sigList$Junction$exprs_newGeneID, sigList$Transcript$exprs_EnsemblGeneID,
	sigList$DER$exprs_geneid[sigDer$exprs_distance == 0]))
allGenesID = allGenesID[!is.na(allGenesID) & !grepl("-", allGenesID)]
length(allGenesID)

convergeMat = matrix(FALSE, nrow = length(allGenesSym), ncol = 5)	
rownames(convergeMat) = allGenesSym
colnames(convergeMat) = c("Gene", "Exon", "Junction", "Transcript", "DER")
convergeMat[which(allGenesSym %in% sigList$Gene$exprs_symbol),"Gene"] = TRUE
convergeMat[which(allGenesSym %in% sigList$Exon$exprs_symbol),"Exon"] = TRUE
convergeMat[which(allGenesSym %in% sigList$Junction$exprs_newGeneSym),"Junction"] = TRUE
convergeMat[which(allGenesSym %in% sigList$Transcript$exprs_ref_gene_id),"Transcript"] = TRUE
convergeMat[which(allGenesSym %in% sigList$DER$exprs_symbol),"DER"] = TRUE
colSums(convergeMat)
table(rowSums(convergeMat))

library(gplots)
pdf("plots/vennDiagram_annotated.pdf")
par(mar=c(1,1,1,1))
venn(as.data.frame(convergeMat))
venn(as.data.frame(convergeMat[,-4])) # no tx
venn(as.data.frame(convergeMat[,-5])) # no der
dev.off()

##### fetal versus control_13plus
pdf("plots/fetal_vs_control_13plus_eqtls.pdf")
par(mar=c(5,6,2,2))
r = quantile(sigGene$beta, c(0.005, 0.995))
plot(beta ~ fetal_slope, data=sigGene,
	xlab="Fetal Estimate", ylab="Adult Estimate",
	main = "Gene-level eQTLs",pch=21,bg="grey",
	ylim = r, xlim=r,cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
abline(0,1,lty=2,col="red", lwd=2)
legend("topleft", paste0("r=",	signif(cor(sigGene$beta,
	sigGene$fetal_slope),3)),	cex=2)
	
r = quantile(sigExon$beta, c(0.005, 0.995))
plot(beta ~ fetal_slope, data=sigExon,
	xlab="Fetal Estimate", ylab="Adult Estimate",
	main = "Exon-level eQTLs",pch=21,bg="grey",
	ylim = r, xlim=r,cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
abline(0,1,lty=2,col="red", lwd=2)
legend("topleft", paste0("r=",	signif(cor(sigExon$beta,
	sigExon$fetal_slope),3)),	cex=2)
	
r = quantile(sigTrans$beta, c(0.005, 0.995))
plot(beta ~ fetal_slope, data=sigTrans,
	xlab="Fetal Estimate", ylab="Adult Estimate",
	main = "Transcript-level eQTLs",pch=21,bg="grey",
	ylim = r, xlim=r,cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
abline(0,1,lty=2,col="red", lwd=2)
legend("topleft", paste0("r=",	signif(cor(sigTrans$beta,
	sigTrans$fetal_slope),3)),	cex=2)
	
r = quantile(sigJxn$beta, c(0.005, 0.995))
plot(beta ~ fetal_slope, data=sigJxn,
	xlab="Fetal Estimate", ylab="Adult Estimate",
	main = "Junction-level eQTLs",pch=21,bg="grey",
	ylim = r, xlim=r,cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
abline(0,1,lty=2,col="red", lwd=2)
legend("topleft", paste0("r=",	signif(cor(sigJxn$beta,
	sigJxn$fetal_slope),3)),	cex=2)
	
dev.off()

##### CAUC vs AA
pdf("plots/cauc_vs_aa.pdf")
par(mar=c(5,6,2,2))
r = quantile(c(sigGene$AA_slope,sigGene$CAUC_slope),
	c(0.005, 0.995))
plot(AA_slope ~ CAUC_slope, data=sigGene,
	xlab="CAUC Estimate", ylab="AA Estimate",
	main = "Gene-level eQTLs",pch=21,bg="grey",
	ylim = r, xlim=r,cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
abline(0,1,lty=2,col="red", lwd=2)
legend("topleft", paste0("r=",	signif(cor(sigGene$AA_slope,
	sigGene$CAUC_slope),3)),	cex=2)
	
r = quantile(c(sigExon$AA_slope,sigExon$CAUC_slope),
	c(0.005, 0.995))
plot(AA_slope ~ CAUC_slope, data=sigExon,
	xlab="CAUC Estimate", ylab="AA Estimate",
	main = "Exon-level eQTLs",pch=21,bg="grey",
	ylim = r, xlim=r,cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
abline(0,1,lty=2,col="red", lwd=2)
legend("topleft", paste0("r=",	signif(cor(sigExon$AA_slope,
	sigExon$CAUC_slope),3)),	cex=2)
	
r = quantile(c(sigTrans$AA_slope,sigTrans$CAUC_slope),
	c(0.005, 0.995))
plot(AA_slope ~ CAUC_slope, data=sigTrans,
	xlab="CAUC Estimate", ylab="AA Estimate",
	main = "Transcript-level eQTLs",pch=21,bg="grey",
	ylim = r, xlim=r,cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
abline(0,1,lty=2,col="red", lwd=2)
legend("topleft", paste0("r=",	signif(cor(sigTrans$AA_slope,
	sigTrans$CAUC_slope),3)),	cex=2)
	
r = quantile(c(sigJxn$AA_slope,sigJxn$CAUC_slope),
	c(0.005, 0.995))
plot(AA_slope ~ CAUC_slope, data=sigJxn,
	xlab="CAUC Estimate", ylab="AA Estimate",
	main = "Junction-level eQTLs",pch=21,bg="grey",
	ylim = r, xlim=r,cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
abline(0,1,lty=2,col="red", lwd=2)
legend("topleft", paste0("r=",	signif(cor(sigJxn$AA_slope,
	sigJxn$CAUC_slope),3)),	cex=2)
	
dev.off()