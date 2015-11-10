setwd("~/Desktop/Super-Enhancers")

#import enhancer list and noodle doc, label columns of noodle doc
enhancer <- read.csv("testEnhancer.csv")

chr7 <- read.table("Noodles/noodles.chr7.txt", quote="", na.strings="", comment.char="",
                   stringsAsFactors=F, header=F,sep="\t",row.names=1)

noodleLabels <- read.table("noodleLabels.csv", header = F,sep=",",
                           quote="", na.strings="", comment.char="", stringsAsFactors=F,
                           nrows=1,row.names=1)
colnames(chr7) <- noodleLabels

#find noodles overlapping with enhancer
noodleOverlap <- chr7[which(chr7$end >= enhancer$start & chr7$start <= enhancer$end), ]

#create plot with both normal & tumor meth ratios
x <- cbind(noodleOverlap$start, noodleOverlap$start)
y <- cbind(noodleOverlap$nor.ratio, noodleOverlap$tmr.ratio)
matplot(x,y, type="l")
legend('topleft',pch='-',col=c('red','black'),legend=c('T','N'))

#find indices of hotspot (max nor.ratio) and determine which samples have 5 or fewer reads
hotspot <- which(max(noodleOverlap$nor.ratio)==noodleOverlap$nor.ratio)
hotspotSamp <- names(which(apply(noodleOverlap[hotspot,grep('700.reads',colnames(chr7))] <= 5,2,all)))
DNAsamp <- names(which(apply(noodleOverlap[hotspot,sub('.700.reads','',hotspotSamp)] == 0,2,all)))

#import clinical data and RNA expression (log)
RNAExp <- read.table("HPVOPRSEMData_091615_log.csv", header = T, sep=",", quote="",
                     na.strings="", comment.char="", stringsAsFactors=F)
RNAExp[2:3] <- list(NULL)
RNAExp[74:75] <- list(NULL)
clinical <- read.table("clinicaldata.csv", header = T,
                       sep=",", quote="", na.strings="", comment.char="", stringsAsFactors=F)

#get expression data for closest gene 
closestGene <- as.character(enhancer$closest.TSS)
closestGene <- paste0(closestGene, "|")
geneExp <- RNAExp[grep(closestGene, RNAExp$genes,fixed=T), ]

#match DNA sample #s to RNA sample #
RNAsamp <- clinical[which(clinical$Tumor.Tissue.DNA.HAND.ID %in% DNAsamp), ]
RNAsamp <- RNAsamp$Tumor.Tissue.RNA.HAND.ID

#separate normal and tumor samples (gene expression data)
normSamp <- geneExp[ ,grep('N',colnames(geneExp))]
tmrSamp <- geneExp[ ,grep('T',colnames(geneExp))]

#choose hypometh samples from N and T groups 
