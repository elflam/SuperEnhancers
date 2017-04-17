###Modificaiton of DefineSamplesAltered.R to remove sample separation by methylation status
#3/24/2017
#Emily Flam

#### support functions 

#serially merge multiple matrices into single matrix 
expmat <- function(x, y){
  df <- merge(x, y, by= "row.names", all.x= F, all.y= F, sort=F)
  rownames(df) <- df$Row.names
  df$Row.names <- NULL
  return(df)
}

###main code
setwd("~/Desktop/Super-Enhancers/Code Skeleton")

#import enhancer list; add rownames to enhancer list including chromosome and genomic coordinates
enhancer <- read.csv("Enhancer Data/2017-03-31FDRDifferentialMethylationEnhancers.csv", header=TRUE, sep=",",quote="",na.strings="", comment.char="", stringsAsFactors = FALSE)
for (i in 1:nrow(enhancer)){
  rownames(enhancer)[i] <- paste(as.character(enhancer[i,]$chr), ":", enhancer[i,]$start, "-", enhancer[i,]$end, sep="")
}

#Load noodles
#noodles.C.Rda contains methylation data
load("/Volumes/Seagate Backup Plus Drive/Noodles.C.Rda")
#huge.report.frame contains coordinate and closest TSS info
load("/Volumes/Seagate Backup Plus Drive/Noodles.huge.report.frame.Rda")

#assign report frame rownames to noodles.C
rownames(noodles.C.methylation) <- rownames(huge.report.frame)

#import RNA expression (log)
#prepare RNAExp for analysis: split the gene list to only include gene names and apply as rownames  
RNAExp <- read.table("ScreenedExpData.csv", header = TRUE, sep=",", quote = "",na.strings="", comment.char="", stringsAsFactors=FALSE, check.names = FALSE)
geneNames <- sapply(strsplit(RNAExp$genes, split="\\|"), function(x){x[[1]]})
rownames(RNAExp) <- geneNames
RNAExp2 <- RNAExp[,-1]

#import clinical data: reference to link RNA sample IDs and meth sample IDs
clinical <- read.table("Enhancer Data/clinicaldata.csv", header = T,sep=",", quote="", na.strings="", comment.char="", stringsAsFactors=F)

#initialize empty list to store information for each enhancer
j=1
enhancerList <- list()

#begin loop #1
for (i in 1:nrow(enhancer)){
  k <- i
  
  #Collect noodles (huge.report.frame and noodles.C.methylation) that overlap with enhancer
  noodleOverlap <- huge.report.frame[which(huge.report.frame$end >= enhancer[i,]$start & huge.report.frame$start <= enhancer[i,]$end & huge.report.frame$chr == enhancer[i,]$chr), ]
  
  noodleMethOverlap <- noodles.C.methylation[which(huge.report.frame$end >= enhancer[i,]$start & huge.report.frame$start <= enhancer[i,]$end & huge.report.frame$chr == enhancer[i,]$chr), ]
  
  #find indices of hotspot where methylation difference between normal and tumor is greatest
  hotspot <- which((abs(noodleOverlap$nor.ratio - noodleOverlap$tmr.ratio) == max(abs(noodleOverlap$nor.ratio - noodleOverlap$tmr.ratio))) & (abs(noodleOverlap$nor.ratio - noodleOverlap$tmr.ratio) > 0))
  
  #Expand enhancers by 1Mbp on each side
  #All genes with TSS within this range will be considered potential target genes
  expStart <- enhancer[i,]$end - 1000000
  expEnd <- enhancer[i,]$start + 1000000
  
  #find noodles (huge.report.frame) overlapping extended enhancer & get list of closest genes and their genomic positions
  extendedOverlap <- huge.report.frame[which(huge.report.frame$end >= expStart & huge.report.frame$start <= expEnd & huge.report.frame$chr == enhancer[i,]$chr), ]
  geneList <- unique(extendedOverlap$closest.TSS)
  
  #get methylation data for extended noodles
  noodleMethExtendedOverlap <- noodles.C.methylation[which(huge.report.frame$end >= expStart & huge.report.frame$start <= expEnd & huge.report.frame$chr == enhancer[i,]$chr),]
  
  #get expression data for closest genes from RNAExp file
  #Genes that are not listed in RNAExp are filtered out
  geneExp <- RNAExp2[0,]
  screenedGenes <- intersect(geneList, geneNames)
  geneExp <- RNAExp2[screenedGenes, ]
  
  #separate normal and tumor samples (gene expression data)
  Nsamp <- geneExp[ ,grep('N',colnames(geneExp))]
  Tsamp <- geneExp[ ,grep('T',colnames(geneExp))]
  
  #change colnames to only display sample HandID
  colnames(Nsamp) <- c(sapply(strsplit(colnames(Nsamp), split="\\."), function(x){x[[2]]}))
  colnames(Tsamp) <- c(sapply(strsplit(colnames(Tsamp), split="\\."), function(x){x[[2]]}))
  
  #create vector of sample IDs with class type
  classVec <- c(colnames(Nsamp), colnames(Tsamp))
  names(classVec) <- classVec
  classVec[colnames(Nsamp)] <- "N"
  classVec[colnames(Tsamp)] <- "T"
  
  #Create color vector for classVec
  classCols <- classVec
  classCols[colnames(Nsamp)] <- "cyan"
  classCols[colnames(Tsamp)] <- "aquamarine"
  
  #turn classVec into a factor 
  classVec <- factor(classVec, levels=c('N','T'))
  
  #Create matrix of potential target genes and their RNA expression values for tumor and normal samples   
  DEMatrix <- Reduce(expmat, list(Nsamp, Tsamp))
  DEMatrix <- as.matrix(DEMatrix[row.names(DEMatrix)%in%row.names(geneExp),])
  
  #Filter out any genes that will not have a biologically detectable difference in expression
  #All genes with max(exp) - min(exp) < 0.5 are dropped
  DEMatrix <- DEMatrix[apply(log2(DEMatrix+1),1,max) - apply(log2(DEMatrix+1),1,min)>0.5, ,drop=FALSE]
  
  #Create and name position vector to include only those genes from noodles that are also in the RNA Exp file
  positions <- extendedOverlap$pos[which(extendedOverlap$closest.TSS %in% row.names(DEMatrix))]
  names(positions) <- extendedOverlap$closest.TSS[which(extendedOverlap$closest.TSS %in% row.names(DEMatrix))]
  positions <- positions[unique(names(positions))]
  
  #Create and name distance vector to include genes from noodles that are also in RNA Exp file
  distance <- sapply(positions, function(x) x-enhancer[i,]$start) 
  names(distance) <- names(positions)
  
  ##filter distances > 1Mbp
  #Removes all information from various vectors if a gene is found to be outside 1Mbp range
  for (i in names(distance)){
    if (distance[i] < -1000000){
      DEMatrix <- DEMatrix[-(which(rownames(DEMatrix) == i)), ,drop=FALSE] 
      positions <- positions[-(which(names(positions)==i))]
      distance <- distance[-(which(names(distance)==i))]
    }
    else if (distance[i] > 1000000){
      DEMatrix <- DEMatrix[-(which(rownames(DEMatrix) == i)), ,drop=FALSE]
      positions <- positions[-(which(names(positions)==i))]
      distance <- distance[-(which(names(distance)==i))]
    }
  }
  
  #add necessary items to jth element of enhancerList and name for enhancer 
  #increase value of j so the next enhancer can fill the next spot in the list
  enhancerList[[j]] <- list(DEMatrix, classVec, as.character(enhancer[k,]$chr),enhancer[k,]$start, enhancer[k,]$end, hotspot, noodleMethOverlap, positions, distance, classCols)  
  j <- j + 1
  
  }

#save all necessary data for next code: DEstatsAltered2.R
save("enhancerList","enhancer", file = "DefineSamplesAltered2.Rda")
  
  
  

