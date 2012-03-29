## BRING IN THE METHYLATION AND EXPRESSION DATA
## LEVEL 3 ONLY
#####

require(synapseClient)

#####
## GRAB THE METHYLATION DATA
#####
methLayer <- loadEntity("167706")
methMat <- methLayer$objects$methMat
methAnn <- methLayer$objects$methAnn

## GET THE NAMES STRAIGHTENED OUT
methSplit <- strsplit(colnames(methMat), "-", fixed=T)
methParticipant <- sapply(methSplit, "[[", 3)
methSample <- substr(sapply(methSplit, "[[", 4), 1, 2)
methVial <- substr(sapply(methSplit, "[[", 4), 3, 3)
methPortion <- substr(sapply(methSplit, "[[", 5), 1, 2)
methAnalyte <- substr(sapply(methSplit, "[[", 5), 3, 3)
methPlate <- sapply(methSplit, "[[", 6)

methNames <- sapply(strsplit(colnames(methMat), "-", fixed=T), function(x){
  blah <- paste(x[1:5], collapse="-")
  blah <- substr(blah, 1, nchar(blah) - 1)
})



#####
## GRAB THE EXPRESSION DATA
#####
exprLayer <- loadEntity("167731")
exprMat <- exprLayer$objects$exprAgilentMat

## GET THE NAMES STRAIGHTENED OUT
exprSplit <- strsplit(colnames(exprMat), "-", fixed=T)
exprParticipant <- sapply(exprSplit, "[[", 3)
exprSample <- substr(sapply(exprSplit, "[[", 4), 1, 2)
exprVial <- substr(sapply(exprSplit, "[[", 4), 3, 3)
exprPortion <- substr(sapply(exprSplit, "[[", 5), 1, 2)
exprAnalyte <- substr(sapply(exprSplit, "[[", 5), 3, 3)
exprPlate <- sapply(exprSplit, "[[", 6)

exprNames <- sapply(strsplit(colnames(exprMat), "-", fixed=T), function(x){
  blah <- paste(x[1:5], collapse="-")
  blah <- substr(blah, 1, nchar(blah) - 1)
  blah
})



table(unique(colnames(meth)) %in% unique(colnames(expr)))





## BRCA1
brca1meth <- rownames(methAnn)[ grep("BRCA1", methAnn$Gene_Symbol) ]
brca1expr <- rownames(exprMat)[ grep("BRCA1", rownames(exprMat)) ]
## BRCA2
brca2meth <- rownames(methAnn)[ grep("BRCA2", methAnn$Gene_Symbol) ]
brca2expr <- rownames(exprMat)[ grep("BRCA2", rownames(exprMat)) ]

