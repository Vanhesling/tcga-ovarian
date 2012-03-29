## TAKE A LOOK AT THE AVAILABLE OVARIAN DATA FROM TCGA
#####

require(synapseClient)
require(affy)
require(snm)

theseNAs <- c("NA", "", " ", "[Not Reported]")


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
methIds <- ovLayers$layer.id[ ovPlatform == "HumanMethylation27" & ovLevel == "Level_3" ]

## LAUNCH INTO GRABBING LEVEL 3 METHYLATION DATA
for( i in methIds ){
  
  tmpEntity <- downloadEntity(i)
  theseFiles <- file.path(tmpEntity$cacheDir, tmpEntity$files[grepl("Methylation", basename(tmpEntity$files))])
  
  for( f in theseFiles ){
    colTmp <- as.character(read.delim(f, nrow=1, header=F, colClasses=c("NULL", "character", "NULL", "NULL", "NULL"), as.is=TRUE))
    tmpDat <- read.delim(f, header=TRUE, colClasses=c("character", "numeric", "NULL", "NULL", "NULL"), as.is=TRUE, skip=1, na.strings=theseNAs)
    rowTmp <- tmpDat[, 1]
    tmpDat <- matrix(cbind(tmpDat[, -1]))
    rownames(tmpDat) <- rowTmp
    colnames(tmpDat) <- colTmp
    if( !exists("methMat") ){
      methMat <- tmpDat
    } else{
      methMat <- cbind(methMat, tmpDat)
    }
  }
  
}

## GRAB THE ANNOTATION INFO FROM ONE FILE
methAnn <- read.delim(f, header=TRUE, as.is=TRUE, skip=1, na.strings=theseNAs)
rownames(methAnn) <- methAnn[, 1]
methAnn <- methAnn[, -c(1:2)]

methLayer <- Layer(list(name="Methylation", type="G", parentId="163905"))
methLayer <- createEntity(methLayer)
methLayer <- addObject(methLayer, methMat)
methLayer <- addObject(methLayer, methAnn)
methLayer <- storeEntity(methLayer)
methLayer
## ID 167706









#####
## EXPRESSION DATA
#####
exprAgilentIds <- ovLayers$layer.id[ grepl("AgilentG4502A", ovPlatform) & ovLevel == "Level_3" ]

myAgilentBatch <- character()
for( i in exprAgilentIds ){
  
  tmpEntity <- downloadEntity(i)
  theseFiles <- file.path(tmpEntity$cacheDir, tmpEntity$files[grepl("U133A", basename(tmpEntity$files))])
  
  for( f in theseFiles ){
    colTmp <- as.character(read.delim(f, nrow=1, header=F, colClasses=c("NULL", "character"), as.is=TRUE))
    tmpDat <- read.delim(f, header=TRUE, colClasses=c("character", "numeric"), as.is=TRUE, skip=1, na.strings=theseNAs)
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







## U133A - TABLED FOR NOW
exprIds <- ovLayers$layer.id[ ovPlatform == "HT_HG-U133A" & ovLevel == "Level_3" ]

myBatch <- character()
## LAUNCH INTO GRABBING LEVEL 1 EXPRESSION DATA (a lot of data)
for( i in exprIds ){
  
  tmpEntity <- downloadEntity(i)
  theseFiles <- file.path(tmpEntity$cacheDir, tmpEntity$files[grepl("U133A", basename(tmpEntity$files))])
  
  for( f in theseFiles ){
    colTmp <- as.character(read.delim(f, nrow=1, header=F, colClasses=c("NULL", "character"), as.is=TRUE))
    tmpDat <- read.delim(f, header=TRUE, colClasses=c("character", "numeric"), as.is=TRUE, skip=1, na.strings=theseNAs)
    rowTmp <- tmpDat[, 1]
    tmpDat <- matrix(cbind(tmpDat[, -1]))
    rownames(tmpDat) <- rowTmp
    colnames(tmpDat) <- colTmp
    if( !exists("exprMat") ){
      exprMat <- tmpDat
    } else{
      exprMat <- cbind(exprMat, tmpDat)
    }
  }
  
  myBatch <- c(myBatch, rep(i, length(theseFiles)))
  
}


## GRAB MAPPING INFORMATION
exprMapEntity <- downloadEntity(ovLayersAll$layer.id[ grepl("mage",ovLayersAll$layer.name) & grepl("U133A", ovLayersAll$layer.name)])
exprMap <- read.delim(file.path(exprMapEntity$cacheDir, exprMapEntity$files[ grepl("sdrf", basename(exprMapEntity$files)) ]), header=T, as.is=T, na.strings=theseNAs)

exprName <- exprMap$Extract.Name
names(exprName) <- exprMap$Array.Data.File
exprName <- exprName[basename(myFiles)]

ab <- ReadAffy(filenames=myFiles, sampleNames=exprName)
rawMat <- exprs(ab)
normMat <- snm(rawMat, adj.var=model.matrix(~factor(myBatch)), int.var=as.data.frame(matrix(1:length(myBatch), ncol=1, nrow=length(myBatch))), rm.adj=TRUE)
exprs(ab) <- normMat

exprFinal <- rma(ab, normalize=FALSE, background=FALSE)


exprLayer <- Layer(list(name="U133A Expression", type="E", parentId="163905"))
exprLayer <- createEntity(exprLayer)
exprLayer <- addObject(exprLayer, exprMat)
exprLayer <- storeEntity(exprLayer)
exprLayer
## ID 






#####
## CLINICAL DATA
#####
clinId <- ovLayersAll$layer.id[ grepl("clinical", ovLayersAll$layer.name, fixed=T) ]
clinLayer <- downloadEntity(clinId)

clin <- lapply(as.list(clinLayer$files), function(x){
  read.delim(file.path(clinLayer$cacheDir, x), header=T, as.is=T, na.strings=theseNAs)
})
names(clin) <- clinLayer$files

lapply(clin, function(x){ grep("uuid", names(x), value=TRUE) })




