## BRING IN THE METHYLATION AND EXPRESSION DATA
## LEVEL 3 ONLY
#####

source("tcgaID.R")

require(synapseClient)
require(IlluminaHumanMethylation27k.db)

options(stringsAsFactors=T)

#####
## GRAB THE METHYLATION DATA
#####
methLayer2 <- loadEntity("168671")
methMat2 <- methLayer2$objects$methMat2
methMat2 <- methMat2[, order(colnames(methMat2))]

methID2 <- tcgaID(id=colnames(methMat2))
methNames2 <- sapply(strsplit(colnames(methMat2), "-", fixed=T), function(x){
  blah <- paste(x[1:4], collapse="-")
  blah <- substr(blah, 1, nchar(blah) - 1)
  blah
})

idm2 <- methID2$Sample == "01"
methMat2 <- methMat2[, idm2]
methID2 <- lapply(methID2, "[", idm2)
methNames2 <- methNames2[idm2]

## THERE IS ONE DUPLICATE - WILL JUST TAKE THE FIRST ONE
dm2 <- !duplicated(methNames2)
methMat2 <- methMat2[, dm2]
methID2 <- lapply(methID2, "[", dm2)
methNames2 <- methNames2[dm2]
colnames(methMat2) <- methNames2


methLayer3 <- loadEntity("168693")
methMat3 <- methLayer3$objects$methMat3
methMat3 <- methMat3[, order(colnames(methMat3))]

methID3 <- tcgaID(id=colnames(methMat3))
methNames3 <- sapply(strsplit(colnames(methMat3), "-", fixed=T), function(x){
  blah <- paste(x[1:4], collapse="-")
  blah <- substr(blah, 1, nchar(blah) - 1)
  blah
})

idm3 <- methID3$Sample == "01"
methMat3 <- methMat3[, idm3]
methID3 <- lapply(methID3, "[", idm3)
methNames3 <- methNames3[idm3]

## THERE IS ONE DUPLICATE - WILL JUST TAKE THE FIRST ONE
dm3 <- !duplicated(methNames3)
methMat3 <- methMat3[, dm3]
methID3 <- lapply(methID3, "[", dm3)
methNames3 <- methNames3[dm3]
colnames(methMat3) <- methNames3

methAnn3 <- methLayer3$objects$methAnn3

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
exprMat <- exprMat[, ide]
exprID <- lapply(exprID, "[", ide)
exprNames <- exprNames[ide]
colnames(exprMat) <- exprNames


#####
## FIND OVERLAP
#####
meth2 <- methMat2
meth3 <- methMat3
expr <- exprMat[, colnames(meth2)]



## GRAB METHYLATION ANNOTATION
methMap <- as.list(IlluminaHumanMethylation27kALIAS2PROBE)

## BRCA1
brca1meth <- meth[ methMap[["BRCA1"]], ]
brca1expr <- expr[grep("BRCA1", rownames(expr)), ]
## BRCA2
brca2meth <- meth[ methMap[["BRCA2"]], ]
brca2expr <- expr[grep("BRCA2", rownames(expr)), ]

