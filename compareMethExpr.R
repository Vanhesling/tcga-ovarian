## BRING IN THE METHYLATION AND EXPRESSION DATA
## LEVEL 3 ONLY
#####

require(synapseClient)

methLayer <- loadEntity("164021")
methMat <- methLayer$objects$methMat
methAnn <- methLayer$objects$methAnn

methAnn[grepl("BRCA", methAnn$Gene_Symbol), ]

## TAKE A LOOK AT SOME METH DATA
boxplot(as.data.frame(methMat), range=0)








exprLayer <- loadEntity("164023")
exprMat <- exprLayer$objects$exprMat

## CHECK FOR OVERLAP
table(colnames(exprMat) %in% colnames(methMat))

