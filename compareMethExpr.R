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
  blah <- paste(x[1:4], collapse="-")
  blah <- substr(blah, 1, nchar(blah) - 1)
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

# methMat3 <- methLayer3$objects$methMat3
# methMat3 <- methMat3[, order(colnames(methMat3))]
# 
# methID3 <- tcgaID(id=colnames(methMat3))
# methNames3 <- sapply(strsplit(colnames(methMat3), "-", fixed=T), function(x){
#   blah <- paste(x[1:4], collapse="-")
#   blah <- substr(blah, 1, nchar(blah) - 1)
#   blah
# })
# 
# idm3 <- methID3$Sample == "01"
# methMat3 <- methMat3[, idm3]
# methID3 <- lapply(methID3, "[", idm3)
# methNames3 <- methNames3[idm3]
# 
# ## THERE IS ONE DUPLICATE - WILL JUST TAKE THE FIRST ONE
# dm3 <- !duplicated(methNames3)
# methMat3 <- methMat3[, dm3]
# methID3 <- lapply(methID3, "[", dm3)
# methNames3 <- methNames3[dm3]
# colnames(methMat3) <- methNames3


# theseBRCA <- rownames(methAnn3)[grep("BRCA", methAnn3$Gene_Symbol)]
# 
# ## GRAB METHYLATION ANNOTATION
# methMap <- as.list(IlluminaHumanMethylation27kALIAS2PROBE)
# 
# geneMap <- as.list(IlluminaHumanMethylation27kSYMBOL)
# geneMap <- unlist(geneMap)




#####
## GRAB THE EXPRESSION DATA
#####
exprLayer <- loadEntity("168673")
exprMat <- exprLayer$objects$exprAgilentMat
exprMat <- exprMat[, order(colnames(exprMat))]

exprID <- tcgaID(id=colnames(exprMat))
exprNames <- sapply(strsplit(colnames(exprMat), "-", fixed=T), function(x){
  blah <- paste(x[1:4], collapse="-")
  blah <- substr(blah, 1, nchar(blah) - 1)
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

## BRCA1
brca1meth <- methAll[rownames(methAnn3)[grep("BRCA1", methAnn3$Gene_Symbol)], ]
brca1expr <- exprAll[grep("BRCA1", rownames(exprAll)), ]
## BRCA2
brca2meth <- methAll[rownames(methAnn3)[grep("BRCA2", methAnn3$Gene_Symbol)], ]
brca2expr <- exprAll[grep("BRCA2", rownames(exprAll)), ]


#####
## CALCULATE 4 DIFFERENT METRICS - STORE IN A LIST FOR BRCA1
#####
brca1Tests <- list()


## TEST 1
## MEAN BETA IN NORMALS
brca1Tests$test1 <- rowMeans(brca1meth[, tumor==0])
names(brca1Tests$test1) <- rownames(brca1meth)

## TEST 2
## DIFFERENCE IN BETA BETWEEN MEAN NORMALS AND 90%ILE OF TUMORS
brca1Tests$test2 <- sapply(as.list(rownames(brca1meth)), function(x){
  quantile(brca1meth[x, tumor==1], probs=0.90, na.rm=T) - mean(brca1meth[x, tumor==0])
})
names(brca1Tests$test2) <- rownames(brca1meth)

## TEST 3
## EXPRESSION FC BETWEEN NORMALS AND TOP 10% TUMORS W/ HIGHEST METHYLATION (BETA)
brca1Tests$test3 <- sapply(as.list(rownames(brca1meth)), function(x){
  tmp <- brca1meth[x, tumor == 1]
  theseTumors <- names(tmp)[ tmp >= quantile(tmp, probs=0.90, na.rm=T) ]
  expr <- c(brca1expr[theseTumors], brca1expr[tumor==0])
  tum <- c(rep(1, length(theseTumors)), rep(0, sum(tumor==0)))
  tmpFit <- lm(expr ~ tum)
  2^(-1*tmpFit$coefficients["tum"])
})
names(brca1Tests$test3) <- rownames(brca1meth)

## TEST 4
## CORRELATION BETWEEN METHYLATION AND EXPRESSION VALUES
brca1Tests$test4 <- sapply(as.list(1:nrow(brca1meth)), function(x){ cor(brca1meth[x, ], brca1expr, method="spearman", use="complete.obs")})
names(brca1Tests$test4) <- rownames(brca1meth)

## RELAXED AND STRINGENT THRESHOLDS FROM TCGA PAPER
brca1Tests$relaxed <- (brca1Tests$test1 < 0.5) + (brca1Tests$test2 > 0.1) + (brca1Tests$test3 > 2) + (brca1Tests$test4 < (-.2))
brca1Tests$stringent <- (brca1Tests$test1 < 0.4) + (brca1Tests$test2 > 0.3) + (brca1Tests$test3 > 3) + (brca1Tests$test4 < (-.3))



