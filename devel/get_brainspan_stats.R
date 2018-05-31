###
library(jaffelab)
library(GenomicRanges)
library(limma)

## load existing stats from LIBD
load("rdas/devStats_controlSamples.rda")

## load brainspan DLPFC data
load("/dcl01/lieber/ajaffe/psychENCODE_Data/BrainSpan/DLPFC/rpkmCounts_brainSpan_DLPFC.rda")

# age groups
pd$ageGroup = cut(pd$Age, c(-1,0,1,10,20,50))

# age spline
fetal = ifelse(pd$Age < 0, 1,0)
birth = pd$Age
birth[birth < 0] = 0 # linear spline
infant = pd$Age - 1
infant[infant < 0] = 0 # linear spline
child = pd$Age - 10
child[child < 0] = 0 # linear spline
teen = pd$Age - 20
teen[teen < 0] = 0 # linear spline

## modeling
mod = model.matrix(~Age + fetal + birth + infant +
	child + teen +  Sex, data=pd)
mod0 = model.matrix(~  Sex, data=pd)

# genomic data
xformData = list(Gene = log2(geneRpkm+1), Exon = log2(exonRpkm+1), 
	Transcript =  log2(tFpkm+1), Junction = log2(jRpkm+1), 
	DER = log2(regionMat+1))
xformData = lapply(xformData, as.matrix)	

# model
statListSpan = lapply(xformData, function(y) {
	cat(".")
	fit = lmFit(y, mod)
	eb = ebayes(fit)
	fit0 = lmFit(y, mod0)
	ff = getF(fit,fit0, y)
	yResid = cleaningY(y, mod, P=7)
	ff$ageCorr = cor(t(yResid), pd$Age) # for directionality switch
	ff$meanRpkm = 2^rowMeans(y)-1
	return(ff)
})

## add fetal effect
modFetal = model.matrix(~ifelse(pd$Age < 0, 1, 0))
statListSpan = mapply(function(y, dat) {
	cat(".")
	fit = lmFit(y, modFetal)
	dat$log2FC_fetal = fit$coef[,2]
	eb = ebayes(fit)
	dat$pval_fetal = eb$p[,2]
	return(dat)
}, xformData, statListSpan,SIMPLIFY=FALSE)

####### compare effects
statList = mapply(function(x,y) {
	cat(".")
	m = match(names(x), rownames(y))
	x$fstat_span = y$fstat[m]
	x$f_pval_span = y$f_pval[m]
	x$ageCorr_span = y$ageCorr[m]
	x$log2FC_fetal_span = y$log2FC_fetal[m]
	x$pval_fetal_span = y$pval_fetal[m]
	x$meanRpkm_span = y$meanRpkm[m]
	return(x)
}, statList, statListSpan, SIMPLIFY=FALSE)

## filter to sig
statListSig = endoapply(statList, 
	function(x) x[which(x$p_bonf < 0.05)])
sapply(statListSig,length)

save(statListSig, file="rdas/signif_dev_stats_plusBrainSpan.rda")