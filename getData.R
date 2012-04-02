## TAKE A LOOK AT THE AVAILABLE OVARIAN DATA FROM TCGA
#####

require(synapseClient)

options(stringsAsFactors=FALSE)
theseNAs <- c("NA", "", " ", "[Not Reported]", "NULL", "null")


tcgaData <- synapseQuery('SELECT * FROM dataset WHERE dataset.repository == "TCGA" AND dataset.parentId == "164427"')
ovId <- tcgaData$dataset.id[ grepl("Ovarian", tcgaData$dataset.name) ]

## GRAB THE OVARIAN LAYERS
ovLayersAll <- synapseQuery(paste('SELECT id, name FROM layer WHERE layer.parentId == "', ovId, '"', sep=""))
ovLayers <- ovLayersAll[ grepl("OV.", ovLayersAll$layer.name, fixed=T), ]

## PARSE THROUGH THE LAYERS TO DETERMINE PLATFORM AND LEVEL OF DATA
ovSplit <- strsplit(ovLayers$layer.name, "OV.", fixed=T)
ovSplit <- sapply(ovSplit, "[[", 2)
ovSplit <- strsplit(ovSplit, ".", fixed=T)

ovPlatform <- sapply(ovSplit, "[[", 1)
ovLevel <- sapply(ovSplit, "[[", 2)


#####
## METHYLATION DATA
#####

## LEVEL 2 DATA
methIds2 <- ovLayers$layer.id[ ovPlatform == "HumanMethylation27" & ovLevel == "Level_2" ]

## LAUNCH INTO GRABBING LEVEL 3 METHYLATION DATA
for( i in methIds2 ){
  
  tmpEntity <- downloadEntity(i)
  theseFiles <- file.path(tmpEntity$cacheDir, tmpEntity$files[grepl("Methylation", basename(tmpEntity$files))])
  
  for( f in theseFiles ){
    colTmp <- as.character(read.delim(f, nrow=1, header=F, colClasses=c("NULL", "character", "NULL", "NULL"), as.is=TRUE))
    tmpDat <- read.delim(f, header=TRUE, colClasses=c("character", "numeric", "NULL", "NULL"), as.is=TRUE, skip=1, na.strings=theseNAs)
    rowTmp <- tmpDat[, 1]
    tmpDat <- matrix(cbind(tmpDat[, -1]))
    rownames(tmpDat) <- rowTmp
    colnames(tmpDat) <- colTmp
    if( !exists("methMat2") ){
      methMat2 <- tmpDat
    } else{
      methMat2 <- cbind(methMat2, tmpDat)
    }
  }
  
}


methLayer2 <- Layer(list(name="Methylation - Level 2", type="G", parentId="163905"))
methLayer2 <- createEntity(methLayer2)
methLayer2 <- addObject(methLayer2, methMat2)
methLayer2 <- storeEntity(methLayer2)
methLayer2
## ID 168671


## LEVEL 3 DATA
methIds3 <- ovLayers$layer.id[ ovPlatform == "HumanMethylation27" & ovLevel == "Level_3" ]

## LAUNCH INTO GRABBING LEVEL 3 METHYLATION DATA
for( i in methIds3 ){
  
  tmpEntity <- downloadEntity(i)
  theseFiles <- file.path(tmpEntity$cacheDir, tmpEntity$files[grepl("Methylation", basename(tmpEntity$files))])
  
  for( f in theseFiles ){
    colTmp <- as.character(read.delim(f, nrow=1, header=F, colClasses=c("NULL", "character", "NULL", "NULL", "NULL"), as.is=TRUE))
    tmpDat <- read.delim(f, header=TRUE, colClasses=c("character", "numeric", "NULL", "NULL", "NULL"), as.is=TRUE, skip=1, na.strings=theseNAs)
    rowTmp <- tmpDat[, 1]
    tmpDat <- matrix(cbind(tmpDat[, -1]))
    rownames(tmpDat) <- rowTmp
    colnames(tmpDat) <- colTmp
    if( !exists("methMat3") ){
      methMat3 <- tmpDat
    } else{
      methMat3 <- cbind(methMat3, tmpDat)
    }
  }
  
}

methAnn3 <- read.delim(f, header=T, as.is=T, skip=1, na.strings=theseNAs)
rowTmp <- methAnn3[, 1]
methAnn3 <- methAnn3[, -c(1,2)]
rownames(methAnn3) <- rowTmp

methLayer3 <- Layer(list(name="Methylation - Level 3", type="G", parentId="163905"))
methLayer3 <- createEntity(methLayer3)
methLayer3 <- addObject(methLayer3, methMat3)
methLayer3 <- addObject(methLayer3, methAnn3)
methLayer3 <- storeEntity(methLayer3)
methLayer3
## ID 168693







#####
## EXPRESSION DATA
#####
exprAgilentIds <- ovLayers$layer.id[ grepl("AgilentG4502A", ovPlatform) & ovLevel == "Level_3" ]

myAgilentBatch <- character()
for( i in exprAgilentIds ){
  
  tmpEntity <- downloadEntity(i)
  theseFiles <- file.path(tmpEntity$cacheDir, tmpEntity$files[grepl("gene", basename(tmpEntity$files))])
  
  for( f in theseFiles ){
    colTmp <- as.character(read.delim(f, nrow=1, header=F, colClasses=c("NULL", "character"), as.is=TRUE))
    tmpDat <- read.delim(f, header=FALSE, colClasses=c("character", "numeric"), as.is=TRUE, skip=2, na.strings=theseNAs)
    rowTmp <- tmpDat[, 1]
    tmpDat <- matrix(cbind(tmpDat[, -1]))
    rownames(tmpDat) <- rowTmp
    colnames(tmpDat) <- colTmp
    if( !exists("exprAgilentMat") ){
      exprAgilentMat <- tmpDat
    } else{
      exprAgilentMat <- cbind(exprAgilentMat, tmpDat)
    }
  }
  
  myAgilentBatch <- c(myAgilentBatch, rep(i, length(theseFiles)))
  
}

exprLayer <- Layer(list(name="Agilent Expression - Level 3", type="E", parentId="163905"))
exprLayer <- createEntity(exprLayer)
exprLayer <- addObject(exprLayer, exprAgilentMat)
exprLayer <- addObject(exprLayer, myAgilentBatch)
exprLayer <- storeEntity(exprLayer)
exprLayer
## ID 168673



#####
## CLINICAL DATA
#####
# clinId <- ovLayersAll$layer.id[ grepl("clinical", ovLayersAll$layer.name, fixed=T) & !grepl("intgen", ovLayersAll$layer.name, fixed=T) ]
# clinL <- downloadEntity(clinId)
# 
# clin <- lapply(as.list(clinL$files), function(x){
#   read.delim(file.path(clinL$cacheDir, x), header=T, as.is=T, na.strings=theseNAs)
# })
# names(clin) <- clinL$files
# 
# clinLayer <- Layer(list(name="Clinical Features", type="C", parentId="163905"))
# clinLayer <- createEntity(clinLayer)
# clinLayer <- addObject(clinLayer, clin)
# clinLayer <- storeEntity(clinLayer)
# clinLayer

#lapply(clin, function(x){ grep("uuid", names(x), value=TRUE) })




