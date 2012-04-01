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
methLayer <- loadEntity("167945")
methMat <- methLayer$objects$methMat
methMat <- methMat[, order(colnames(methMat))]

methID <- tcgaID(id=colnames(methMat))
methNames <- sapply(strsplit(colnames(methMat), "-", fixed=T), function(x){
  blah <- paste(x[1:4], collapse="-")
  blah <- substr(blah, 1, nchar(blah) - 1)
  blah
})

idm <- methID$Sample == "01"
methMat <- methMat[, idm]
methID <- lapply(methID, "[", idm)
methNames <- methNames[idm]

## THERE IS ONE DUPLICATE - WILL JUST TAKE THE FIRST ONE
dm <- !duplicated(methNames)
methMat <- methMat[, dm]
methID <- lapply(methID, "[", dm)
methNames <- methNames[dm]
colnames(methMat) <- methNames

#####
## GRAB THE EXPRESSION DATA
#####
exprLayer <- loadEntity("167731")
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
meth <- methMat
expr <- exprMat[, colnames(meth)]

## GRAB METHYLATION ANNOTATION
methMap <- as.list(IlluminaHumanMethylation27kALIAS2PROBE)

## BRCA1
brca1meth <- meth[ methMap[["BRCA1"]], ]
brca1expr <- expr[grep("BRCA1", rownames(expr)), ]
## BRCA2
brca2meth <- meth[ methMap[["BRCA2"]], ]
brca2expr <- expr[grep("BRCA2", rownames(expr)), ]

