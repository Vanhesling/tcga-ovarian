## BRING IN THE METHYLATION AND EXPRESSION DATA
## LEVEL 3 ONLY
#####

require(synapseClient)

methLayer <- loadEntity("")
methMat <- methLayer$objects$methMat
methAnn <- methLayer$objects$methAnn


## BRCA1
brca1meth <- rownames(methAnn)[ grep("BRCA1", methAnn$Gene_Symbol) ]
## BRCA2
brca2meth <- rownames(methAnn)[ grep("BRCA2", methAnn$Gene_Symbol) ]


## TAKE A LOOK AT SOME METH DATA
boxplot(as.data.frame(methMat), range=0)








exprLayer <- loadEntity("164023")
exprMat <- exprLayer$objects$exprMat

## CHECK FOR OVERLAP
table(colnames(exprMat) %in% colnames(methMat))

