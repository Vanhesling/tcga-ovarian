## BRING IN THE METHYLATION AND EXPRESSION DATA
## ONLY PUBLICALLY AVAILABLE
#####

source("tcgaID.R")

require(synapseClient)
#require(IlluminaHumanMethylation27k.db)

options(stringsAsFactors=T)

#####
## GRAB THE METHYLATION DATA
#####

## USE THE LEVEL 2 DATA - BUT PULL IN THE ANNOTATION INFO FROM LEVEL 3
methLayer <- loadEntity("168671")
methMat <- methLayer$objects$methMat2
methMat <- methMat[, order(colnames(methMat))]

methID <- tcgaID(id=colnames(methMat))
methNames <- sapply(strsplit(colnames(methMat), "-", fixed=T), function(x){
  blah <- paste(x[1:3], collapse="-")
  blah
})

idm <- methID$Sample == "01"
methTumor <- methMat[, idm]
methTumorID <- lapply(methID, "[", idm)
methTumorNames <- methNames[idm]

idmn <- methID$Sample %in% as.character(10:19) & methID$TSS == "01" & methID$Portion == "01"
methNormal <- methMat[, idmn]
methNormalID <- lapply(methID, "[", idmn)
methNormalNames <- methNames[idmn]
colnames(methNormal) <- methNormalNames

## THERE IS ONE DUPLICATE - WILL JUST TAKE THE FIRST ONE
dm <- !duplicated(methTumorNames)
methTumor <- methTumor[, dm]
methTumorID <- lapply(methTumorID, "[", dm)
methTumorNames <- methTumorNames[dm]
colnames(methTumor) <- methTumorNames


## JUST GRAB THE ANNOTATION FROM THE LEVEL 3 DATA
methLayer3 <- loadEntity("168693")
methAnn3 <- methLayer3$objects$methAnn3



#####
## GRAB THE EXPRESSION DATA
#####
exprLayer <- loadEntity("168673")
exprMat <- exprLayer$objects$exprAgilentMat
exprMat <- exprMat[, order(colnames(exprMat))]

exprID <- tcgaID(id=colnames(exprMat))
exprNames <- sapply(strsplit(colnames(exprMat), "-", fixed=T), function(x){
  blah <- paste(x[1:3], collapse="-")
  blah
})

ide <- exprID$Sample == "01"
exprTumor <- exprMat[, ide]
exprTumorID <- lapply(exprID, "[", ide)
exprTumorNames <- exprNames[ide]
colnames(exprTumor) <- exprTumorNames

iden <- exprID$Sample %in% as.character(10:19) & exprID$Portion == "01"
exprNormal <- exprMat[, iden]
exprNormalID <- lapply(exprID, "[", iden)
exprNormalNames <- exprNames[iden]
colnames(exprNormal) <- exprNormalNames

#####
## FIND OVERLAP
#####
exprTumor <- exprTumor[, colnames(methTumor)]
exprNormal <- exprNormal[, colnames(methNormal)]

methAll <- cbind(methTumor, methNormal)
exprAll <- cbind(exprTumor, exprNormal)
tumor <- c(rep(1, ncol(methTumor)), rep(0, ncol(methNormal)))






#####
## CALCULATE 4 DIFFERENT METRICS
#####
cpgTests <- function(geneSymbol){
  tmpMeth <- methAll[rownames(methAnn3)[grep(geneSymbol, methAnn3$Gene_Symbol)], ]
  tmpExpr <- exprAll[grep(geneSymbol, rownames(exprAll)), ]
  
  res <- list()
  
  ## TEST 1
  ## MEAN BETA IN NORMALS
  res$test1 <- rowMeans(tmpMeth[, tumor==0])
  names(res$test1) <- rownames(tmpMeth)
  
  ## TEST 2
  ## DIFFERENCE IN BETA BETWEEN MEAN NORMALS AND 90%ILE OF TUMORS
  res$test2 <- sapply(as.list(rownames(tmpMeth)), function(x){
    quantile(tmpMeth[x, tumor==1], probs=0.90, na.rm=T) - mean(tmpMeth[x, tumor==0])
  })
  names(res$test2) <- rownames(tmpMeth)
  
  ## TEST 3
  ## EXPRESSION FC BETWEEN NORMALS AND TOP 10% TUMORS W/ HIGHEST METHYLATION (BETA)
  res$test3 <- sapply(as.list(rownames(tmpMeth)), function(x){
    tmp <- tmpMeth[x, tumor == 1]
    theseTumors <- names(tmp)[ tmp >= quantile(tmp, probs=0.90, na.rm=T) ]
    expr <- c(tmpExpr[theseTumors], tmpExpr[tumor==0])
    tum <- c(rep(1, length(theseTumors)), rep(0, sum(tumor==0)))
    tmpFit <- lm(expr ~ tum)
    2^(-1*tmpFit$coefficients["tum"])
  })
  names(res$test3) <- rownames(tmpMeth)
  
  ## TEST 4
  ## CORRELATION BETWEEN METHYLATION AND EXPRESSION VALUES
  res$test4 <- sapply(as.list(1:nrow(tmpMeth)), function(x){ cor(tmpMeth[x, ], tmpExpr, method="spearman", use="complete.obs")})
  names(res$test4) <- rownames(tmpMeth)
  
  ## RELAXED AND STRINGENT THRESHOLDS FROM TCGA PAPER
  res$relaxed <- (res$test1 < 0.5) + (res$test2 > 0.1) + (res$test3 > 2) + (res$test4 < (-.2))
  res$stringent <- (res$test1 < 0.4) + (res$test2 > 0.3) + (res$test3 > 3) + (res$test4 < (-.3))
  
  return(res)
}

## TRY SLIDING THRESHOLD FOR TOP 10 PERCENT



myRes <- list()
myRes[["BRCA1"]] <- cpgTests("BRCA1")
myRes[["BRCA2"]] <- cpgTests("BRCA2")
myRes[["RAD51C"]] <- cpgTests("RAD51C")
myRes[["RAD51L3"]] <- cpgTests("RAD51L3")

#####
## EPIGENETIC SILENCING -- AS PER TCGA PAPER -- SLIGHTLY DIFFERENT DATA STANDARDIZATION
#####
silenced <- function(cpg){
  tmpMeth <- as.numeric(methTumor[cpg, ])
  tmpExpr <- as.numeric(exprTumor[methAnn3$Gene_Symbol[which(rownames(methAnn3) == cpg)], ])
  
  tmpMeth <- (tmpMeth - mean(tmpMeth, na.rm=T))/sd(tmpMeth, na.rm=T)
  tmpExpr <- (tmpExpr - mean(tmpExpr, na.rm=T))/sd(tmpExpr, na.rm=T)
  
  myMat <- cbind(tmpMeth, tmpExpr)
  rownames(myMat) <- colnames(methTumor)
  myMat <- na.exclude(myMat)
  
  myK <- kmeans(myMat, centers=2)
  silenced <- ifelse(abs(diff(myK$centers[1, ])) > abs(diff(myK$centers[2, ])), 1, 2)
  silenced <- ifelse(myK$cluster == silenced, 1, 0)
  silenced
}

#####
## TEST FOR BRCA1 SILENCING
#####

brca1CpGs <- names(myRes$BRCA1$test4[myRes$BRCA1$test4 < (-0.2)])
# > brca1CpGs
# [1] "cg04658354" "cg08993267" "cg19088651" "cg19531713"
## THESE MATCH THE PUBLISHED 4 CPGS

brca1silenced <- lapply(as.list(brca1CpGs), silenced)
brca1silenced <- unlist(brca1silenced)
cts <- table(names(brca1silenced), brca1silenced)

## PERCENT OF CPGS SHOWING SILENCING
brca1silenced <- as.numeric(cts[, "1"]) / rowSums(cts)

## LIST OF SAMPLES HAVING GREATER THAN 50% SILENCING
silencedSamples <- names(brca1silenced)[ brca1silenced > 0 ]
brca1silenced <- brca1silenced[ brca1silenced > 0 ]




write.table(cbind(patient.id=silencedSamples, gene.silenced=rep("BRCA1", length(silencedSamples)), pct.cpgs=brca1silenced), file="~/epigeneticSilenced.txt", sep="\t", quote=F, row.names=F, col.names=T)
