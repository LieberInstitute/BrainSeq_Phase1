###
###
source("../eqtl_functions.R")

library(GenomicRanges)
library(clusterProfiler)

##### summary statistics for all features ####
load("rdas/devStats_controlSamples.rda")

#### and those with a switch
load("rdas/isoform_switch_devel_byFeature.rda")

####################################
##### gwas enrichment for switch ###
geneMapGR = statList$Gene

## read in GWAS
pgcFinal = read.delim("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/tables/pgc2_loci.txt", as.is=TRUE)
pgcFinal$chr = ss(pgcFinal$Position..hg19.,":")
tmp = ss(pgcFinal$Position..hg19.,":",2)
pgcFinal$start = as.numeric(ss(tmp, "-"))
pgcFinal$end = as.numeric(ss(tmp, "-",2))
gr108 = GRanges(pgcFinal$chr, IRanges(pgcFinal$start,pgcFinal$end))

## compared to expressed background
geneSwitchList = lapply(switchList, rownames)
geneBgList = lapply(statList, function(x) 
	unique(x$EnsemblGeneID[!is.na(x$p_bonf)]))
geneBgList = lapply(geneBgList, function(x)
	x[!is.na(x) & !grepl("-", x)])
elementLengths(geneBgList)

## get expressed genes in each PGC region
pgcExprsList = lapply(statList, function(x) {
	x = x[!is.na(x$p_bonf)]
	oo1 = findOverlaps(gr108, x)
	ens = x$EnsemblGeneID[subjectHits(oo1)]
	unique(ens[!is.na(ens) & !grepl("-",ens)])
})

#### enrichment test
matPgcList = mapply(function(g,gInPgc,bg) {
	g = g[!is.na(g) & !grepl("-", g)]
	m = matrix(c(sum(gInPgc %in% g),
		sum(! gInPgc %in% g), 
		sum(! g %in% gInPgc),
		length(bg) - length(union(g, gInPgc))),
			byrow=2, nc = 2)
	rownames(m) = c("inReg", "notReg")
	colnames(m) = c("isoSw", "noIso")
	m
}, geneSwitchList, pgcExprsList[-1], geneBgList[-1],SIMPLIFY=FALSE)
matPgcList

## p-value and ORs for association
sapply(matPgcList, function(x) chisq.test(x)$p.value)
sapply(matPgcList, function(x) x[1,1]*x[2,2]/x[1,2]/x[2,1])

####### 
## dev

pgcDevList = lapply(statList, function(x) {
	x = x[which(x$p_bonf < 0.05)]
	oo1 = findOverlaps(gr108, x)
	ens = x$EnsemblGeneID[subjectHits(oo1)]
	unique(ens[!is.na(ens) & !grepl("-",ens)])
})

geneBgDevList = lapply(statList, function(x)
	unique(x$EnsemblGeneID[which(x$p_bonf < 0.05)]))
geneBgDevList = lapply(geneBgDevList, function(x) x[!is.na(x) & !grepl("-", x)])
#### enrichment test
matPgcList_dev = mapply(function(g,gInPgc,bg) {
	g = g[!is.na(g) & !grepl("-", g)]
	m = matrix(c(sum(gInPgc %in% g),
		sum(! gInPgc %in% g), 
		sum(! g %in% gInPgc),
		length(bg) - length(union(g, gInPgc))),
			byrow=2, nc = 2)
	rownames(m) = c("inReg", "notReg")
	colnames(m) = c("isoSw", "noIso")
	m
}, geneSwitchList, pgcDevList[-1], geneBgDevList[-1],SIMPLIFY=FALSE)
matPgcList_dev

## p-value and ORs for association
sapply(matPgcList_dev, function(x) chisq.test(x)$p.value)
sapply(matPgcList_dev, function(x) x[1,1]*x[2,2]/x[1,2]/x[2,1])

####### check expression 
statList2 = mapply(function(x,y) {
	cat(".")
	x$hasSwitch = x$EnsemblGeneID %in% rownames(y)
	x$inPgc = countOverlaps(x, gr108) > 0
	return(x)
}, statList[-1], switchList)

statList2 = endoapply(statList2, function(x) x[!is.na(x$p_bonf)])

sapply(statList2, function(x) {
	tapply(x$meanControl, paste0(x$hasSwitch,":", x$inPgc), mean)
})

adjOr = t(sapply(statList2, function(x) {
	summary(glm( x$inPgc ~ x$hasSwitch + 
		x$meanControl,family="binomial"))$coef[2,]
}))

#############################
# compare to null regions ###
#############################
xx = load("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/nullRegions_pgc2_B10000.rda")

## null by type
geneMap = statList$Gene
pgcNullList_all = lapply(statList, function(x) {
	cat(".")
	x = x[!is.na(x$p_bonf)]
	map = geneMap[names(geneMap) %in% x$EnsemblGeneID]
	oo1 = findOverlaps(nullRegions, map)
	CharacterList(split(names(map)[subjectHits(oo1)],
		nullRegions$permutation[queryHits(oo1)]))
})

### enrichment test compared to expressed
matNull_all = mapply(function(g,gInPgc,bg) {
	g = g[!is.na(g) & !grepl("-", g)]
	topleft = sum(gInPgc %in% g)
	topright = sum(! gInPgc %in% g)
	bottomleft = length(g) - topleft
	bottomright = length(bg) - 
		topleft - topright - bottomleft
	topleft/bottomleft/topright*bottomright
}, geneSwitchList, pgcNullList_all[-1], geneBgList[-1])

obsOR = sapply(matPgcList, function(x) x[1,1]*x[2,2]/x[1,2]/x[2,1])
obsORMat = matrix(obsOR, nc = ncol(matNull_all),
	nr = nrow(matNull_all),byrow=TRUE)
empP = colMeans(matNull_all > obsORMat)
empP

pdf("plots/nullDistributions_exprsOverlap_PGC2.pdf",w=5,h=4)
for(i in 1:ncol(matNull_all)) {
	hist(matNull_all[,i], main = colnames(matNull_all)[i],
		col="grey",xlab="OR",cex.axis=2.2)
	abline(v=obsOR[i], col="red",lwd=3)
	legend("topright", paste0("p=",signif(empP,3))[i],cex=2,bty="n")
}
dev.off()


#### dev
pgcNullList_dev = lapply(statList, function(x) {
	cat(".")
	x = x[which(x$p_bonf < 0.05)]
	map = geneMap[names(geneMap) %in% x$EnsemblGeneID]
	oo1 = findOverlaps(nullRegions, map)
	CharacterList(split(names(map)[subjectHits(oo1)],
		nullRegions$permutation[queryHits(oo1)]))
})

### enrichment test compared to devel regulated
matNull_dev = mapply(function(g,gInPgc,bg) {
	g = g[!is.na(g) & !grepl("-", g)]
	topleft = sum(gInPgc %in% g)
	topright = sum(! gInPgc %in% g)
	bottomleft = length(g) - topleft
	bottomright = length(bg) - 
		topleft - topright - bottomleft
	topleft/bottomleft/topright*bottomright
}, geneSwitchList, pgcNullList_dev[-1], geneBgDevList[-1])

obsOR_dev = sapply(matPgcList_dev, function(x) x[1,1]*x[2,2]/x[1,2]/x[2,1])
obsORMat_dev = matrix(obsOR_dev, nc = ncol(matNull_dev),
	nr = nrow(matNull_dev),byrow=TRUE)
colMeans(matNull_dev > obsORMat_dev)

########################
## what are the genes? #
devSwitchGenesPgc = lapply(geneSwitchList, 
	function(g) {
		g = g[!is.na(g) & !grepl("-", g)]
		gg = gInPgc2[gInPgc2 %in% g]
		geneMap[gg,"Symbol"]
})
uGenePgcSwitch = unique(unlist(devSwitchGenesPgc))
pgcMatSwitch = sapply(devSwitchGenesPgc, function(x) uGenePgcSwitch %in% x)
rownames(pgcMatSwitch) = uGenePgcSwitch
pgcMatSwitch = pgcMatSwitch[order(rowSums(pgcMatSwitch),decreasing=TRUE),]
pgcMatSwitch = pgcMatSwitch[rownames(pgcMatSwitch) != "",]
write.csv(pgcMatSwitch, file = "tables/suppTable7_develSwitch_pgcHits.csv")

######
### compare to meQTL genes
meqtls = read.delim("/users/ajaffe/Lieber/Projects/450k/ECD2014/tables/jaffe_suppTable11_pgc2_meQTLs.txt.gz",
	as.is=TRUE)
meqtls = meqtls[meqtls$pvalue < 5e-8,]
meqtlGenes = CharacterList(split(meqtls$meqtlNearestGene,
	meqtls$hitIndex))
meqtlGenes = unlist(lapply(meqtlGenes, unique))
sum(rownames(pgcMatSwitch)%in% meqtlGenes)

#### enrichment
meqtlEnrMat = matrix(FALSE, nc = 2, nr = length(symInPgc),
	dimnames = list(symInPgc, c("isSwitch", "isMeqtl")))
meqtlEnrMat[,"isSwitch"] = symInPgc %in% rownames(pgcMatSwitch)
meqtlEnrMat[,"isMeqtl"] = symInPgc %in% meqtlGenes
tab = table(meqtlEnrMat[,1], meqtlEnrMat[,2],
	dnn = c("isSwitch", "isMeqtl"))
chisq.test(tab) # p=1.96e-8
tab[1,1]*tab[2,2]/tab[1,2]/tab[2,1]

## pgc specific enrichment
entrezPgc = geneMap[gInPgc2, "EntrezID"]
entrezPgc = entrezPgc[!is.na(entrezPgc)]
goListSwitchPgc = lapply(geneSwitchList, function(g) {
	cat(".")
	g_eid = geneMap[gInPgc2[which(gInPgc2 %in% g)],"EntrezID"]
	g_eid = g_eid[!is.na(g_eid)]
	goParams <- new("GOHyperGParams", geneIds = unique(g_eid),
				universeGeneIds = unique(entrezPgc),	
				annotation = "org.Hs.eg.db",
				ontology = "BP", pvalueCutoff = 1, conditional = FALSE,
				testDirection="over")
	ht=hyperGTest(goParams)
	summary(ht)
})


################################
#### compared to non-SZ GWAS ###

xx=load("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/gwasResults_lifted.rda")

## get genes in each region, all expressed first
gwasList = split(gwasLift, gwasLift$Dx)
ooOther_all = findOverlaps(gwasList, geneMapGR_all)
gOther_all = split(names(geneMapGR_all)[subjectHits(ooOther_all)], 
	queryHits(ooOther_all))
names(gOther_all) = names(gwasList)	

matOtherList_all = unlist(lapply(geneSwitchList, function(g) {
	g = g[!is.na(g) & !grepl("-", g)]
	lapply(gOther_all, function(gInOther) {
		matrix(c(sum(gInOther %in% g),
			sum(! gInOther %in% g), 
			sum(! g %in% gInOther),
			length(geneMapGR_all) - length(union(g, gInOther))),
				byrow=2, nc = 2)
	})
}),recur=FALSE)

sapply(matOtherList_all, function(x) chisq.test(x)$p.value)
sapply(matOtherList_all, function(x) x[1,1]*x[2,2]/x[1,2]/x[2,1])


ooOther_dev = findOverlaps(gwasList, geneMapGR_dev)
gOther_dev = split(names(geneMapGR_dev)[subjectHits(ooOther_dev)], 
	queryHits(ooOther_dev))
names(gOther_dev) = names(gwasList)	

matOtherList_dev = unlist(lapply(geneSwitchList, function(g) {
	g = g[!is.na(g) & !grepl("-", g)]
	lapply(gOther_dev, function(gInOther) {
		matrix(c(sum(gInOther %in% g),
			sum(! gInOther %in% g), 
			sum(! g %in% gInOther),
			length(geneMapGR_dev) - length(union(g, gInOther))),
				byrow=2, nc = 2)
	})
}),recur=FALSE)

sapply(matOtherList_dev, function(x) chisq.test(x)$p.value)

sapply(matOtherList_dev, function(x) x[1,1]*x[2,2]/x[1,2]/x[2,1])


