###3/24/2017
###Emily Flam
###Same as DEStatsAltered, but only divides samples into T vs N, no meth status

#required packages
library("edgeR")
library(limma)
library('simpleaffy')

### code
setwd("~/Desktop/Super-Enhancers/Code Skeleton")
load("DefineSamplesAltered2.Rda")

#import RNA expression (log)
#prepare RNAExp for analysis: split the gene list to only include gene names and apply as rownames
RNAExp <- read.table("ScreenedExpData.csv", header = T, sep=",", quote="", na.strings="", comment.char="", stringsAsFactors=F)
geneNames <- sapply(strsplit(RNAExp$genes, split="\\|"), function(x){x[[1]]})
rownames(RNAExp) <- geneNames
RNAExp2 <- RNAExp[,-1]

#import clinical data: reference to link RNA sample IDs and meth sample IDs
clinical <- read.table("Enhancer Data/clinicaldata.csv", header = T, sep=",", quote="", na.strings="", comment.char="", stringsAsFactors=F)

#initialize vector that will be used to order the large enhancer list, putting enhncers with strongest differential expression at the top
#will be filled by manually increase value of j as each enhancer is completed
j<-1
orderVec <- vector()
  
#Create DGElist with all genes, subset RNAExp2 to only include T and N samples
DEMatrix_DGE <- DGEList(counts=RNAExp2[,grep("T|N",colnames(RNAExp2))])
DEMatrix_DGE <- calcNormFactors(DEMatrix_DGE)

for (i in 1:length(enhancerList)){
  
  k<-i

  ##continue differential methylation
  #create factor of T vs N phenotype of each sample 
  pheno <- factor(enhancerList[[i]][[2]][colnames(enhancerList[[i]][[1]])])
  
  #create model matrix from this factor: colnames are "T" and "N"
  mm <- model.matrix(~pheno)
  colnames(mm) <- levels(pheno)
  
  #fit data with model matrix and eBayes
  DEMatrix_VOOM <- voom(DEMatrix_DGE, mm)
  fit <- lmFit(DEMatrix_VOOM,mm)
  DEMatrix.fit <- eBayes(fit) 
  
  #Create contrast matrix, using T and N as opposing groups
  contrast.matrix = makeContrasts(contrasts=c("N - T"), levels=pheno)
  
  #refit data with contrast matrix and eBayes
  fit2= contrasts.fit(DEMatrix.fit, contrast.matrix) 
  fit2= eBayes(fit2) 

  DEgenes <- row.names(topTable(fit2[rownames(enhancerList[[i]][[1]]),],p.value = 0.05, number=Inf,lfc=1))
    
  #prevent 1-row matrices
  if (nrow(enhancerList[[i]][[1]]) == 1){
    name <- rownames(enhancerList[[i]][[1]])
    enhancerList[[i]][[1]] <- rbind(enhancerList[[i]][[1]][1,], enhancerList[[i]][[1]][1,], enhancerList[[i]][[1]][1,])
    rownames(enhancerList[[i]][[1]]) <- c(name, paste(name, "*" ), paste(name,"**"))
    enhancerList[[i]][[8]][2:3] <- c(enhancerList[[i]][[8]][1],enhancerList[[i]][[8]][1])
    enhancerList[[i]][[9]][2:3] <- c(enhancerList[[i]][[9]][1],enhancerList[[i]][[9]][1])
    
   #vector with DE determination
    DE <- row.names(enhancerList[[i]][[1]])
    names(DE) <- row.names(enhancerList[[i]][[1]])
    DE[] <- c("no")
    DE[DEgenes] <- "yes"
    
    #create table with genes, distance to enhancer, location, p-value, tstat
    geneTable <- data.frame(position=enhancerList[[i]][[8]], distance=enhancerList[[i]][[9]], LFC=rep(round(fit2[name,]$coefficients[,1],2),3), p.value = rep(signif(fit2[name,]$p.value[,1],2),3), adj.p.value = rep(signif(p.adjust(fit2[name,]$p.value[,1], method='BH'),4),3), tstat = rep(signif(fit2[name,]$t[,1],2),3), DE = DE)
    
    #for use in LocalMeth.R
    enhancerList[[i]][[20]] <- "X"
    
  }else{
  
    #vector with DE determination
    DE <- row.names(enhancerList[[i]][[1]])
    names(DE) <- row.names(enhancerList[[i]][[1]])
    DE[] <- c("no")
    DE[DEgenes] <- "yes"
  
    #create table with genes, distance to enhancer, location, p-value, tstat
    geneTable <- data.frame(position=enhancerList[[i]][[8]], distance=enhancerList[[i]][[9]], LFC=round(fit2[rownames(enhancerList[[i]][[1]]),]$coefficients[,1],2), p.value = signif(fit2[rownames(enhancerList[[i]][[1]]),]$p.value[,1],2), adj.p.value =signif(p.adjust(fit2[rownames(enhancerList[[i]][[1]]),]$p.value[,1], method='BH'),4), tstat = signif(fit2[rownames(enhancerList[[i]][[1]]),]$t[,1],2), DE = DE)
    
    #for use in LocalMeth.R
    enhancerList[[i]][[20]] <- "O"
    
  }
  
  geneTable <- geneTable[order(geneTable$p.value),]
  
  #sort ExpMatrix to match geneTable order
  enhancerList[[i]][[1]] <- enhancerList[[i]][[1]][rownames(geneTable),]
  
  #add fit2 and DEgenes to enhancer enhancerList list
  enhancerList[[i]][[11]] <- fit2
  enhancerList[[i]][[12]] <- geneTable
  
  
  if (nrow(geneTable) < 5){
    ap <- nrow(geneTable)
  }else{
    ap <- 5
  }
  
  pmean <- mean(geneTable$p.value[1:ap])  
  orderVec[j] <- pmean
  names(orderVec)[j] <- rownames(enhancer)[k]
    
  j <- j+1
}


enhancerList <- enhancerList[order(orderVec)]
sortedEnhancer <- enhancer[names(orderVec[order(orderVec)]),]

save("enhancerList","sortedEnhancer", file = "DEstatsAltered2.Rda" )
