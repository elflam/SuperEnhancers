###3/17/2017
###Emily Flam
###Same as LocalMethKendall, but using only RNAExp columns for T or N samples

# i -> digits; items in enhancerList
# a -> character; genes/rownames of enhancerList[[i]][[1]]
# q -> digits; fill unmethSampsD (manually increase)
#       reset and use to fill newMeth (manually increase)
# m -> digits; fill avgVec, LM, and LM2 (manually increase)
# n -> digits; loop through columns of Poodles to fill meth
# z -> sample #s; loops through members of unmethSampsD2, used to calculate avg for and name newMeth items

setwd("~/Desktop/Super-Enhancers/Code Skeleton")
load("DEstatsAltered2.Rda")

library(Kendall)

#Load noodles
#noodles.C.Rda contains methylation data
load("/Volumes/Seagate Backup Plus Drive/Noodles.C.Rda")
#huge.report.frame contains coordinate and closest TSS info
load("/Volumes/Seagate Backup Plus Drive/Noodles.huge.report.frame.Rda")

#assign report frame rownames to noodles.C
rownames(noodles.C.methylation) <- rownames(huge.report.frame)

#import RNA expression (log)
#prepare RNAExp for analysis  
RNAExp <- read.table("ScreenedExpData.csv", header = T, sep=",", quote="", na.strings="", comment.char="", stringsAsFactors=F)
geneNames <- sapply(strsplit(RNAExp$genes, split="\\|"), function(x){x[[1]]})
rownames(RNAExp) <- geneNames
RNAExp2 <- RNAExp[,-1]
RNAExp2 <- RNAExp2[,grep("T|N",colnames(RNAExp2))]

clinical <- read.table("Enhancer Data/clinicaldata.csv", header = T, sep=",", quote="", na.strings="", comment.char="", stringsAsFactors=F)

#separate normal and tumor samples (gene expression data)
Nsamp <- c(sapply(strsplit(grep("N",colnames(RNAExp), value=TRUE), split="\\."), function(x){x[[2]]}))
Tsamp <- c(sapply(strsplit(grep("T",colnames(RNAExp), value=TRUE), split="\\."), function(x){x[[2]]}))

#Convert phenotype vectors from RNA IDs to meth IDs
Normals <- clinical$Tumor.Tissue.DNA.HAND.ID[which(clinical$Tumor.Tissue.RNA.HAND.ID %in% Nsamp)]
Tumors <- clinical$Tumor.Tissue.DNA.HAND.ID[which(clinical$Tumor.Tissue.RNA.HAND.ID %in% Tsamp)]
TNsamps <- c(Normals,Tumors)

for (i in 1:length(enhancerList)){
  
  m <- 1
  
  #initialize vector for collection of p.vals local promoter meth vs gene exp and average meth for promoter 
  LM <- numeric()
  LMP <- numeric()
  LM2 <- numeric()
  LM2P <- numeric()
  avgMeth <- numeric()
  
  for (a in rownames(enhancerList[[i]][[1]])){
    
    #Create Exp vector for gene
    Exp <- as.numeric(RNAExp2[which(rownames(RNAExp2) == a),])
    names(Exp) <- c(sapply(strsplit(colnames(RNAExp2), split = "\\."), function(x){x[[2]]}))
    Exp <- Exp[sort(names(Exp))]
    
    #get promoter region
    TSS <- enhancerList[[i]][[8]][a] 
    promoterStart <- TSS-1500
    promoterEnd <- TSS + 500
    
    #noodles that span the promoter region
    Poodles <- noodles.C.methylation[which(huge.report.frame$end >= promoterStart & huge.report.frame$start <= promoterEnd & huge.report.frame$chr == enhancerList[[i]][[3]]), which(colnames(noodles.C.methylation) %in% TNsamps)] 
    
    #Create vector of samples with both meth and RNA data
    Allsamps <- colnames(Poodles)
    AllRNAsamps <-clinical$Tumor.Tissue.RNA.HAND.ID[which(clinical$Tumor.Tissue.DNA.HAND.ID %in% Allsamps)]
    
    #Get binary of all noodles per sample
    #Create subset of meth vector for only unmeth samples (based on average of noodles) for later Wilcoxon test
    meth <- numeric()
    unmethSampsD <- character()
    q<-1
    
    for (n in 1:ncol(Poodles)){
      meth[n] <- length(which(Poodles[,n]>0))/length(Poodles[,n])
      names(meth)[n] <- colnames(Poodles)[n]
      if (meth[n] == 0 ){
        unmethSampsD[q] <- colnames(Poodles)[n]
        q <- q + 1 
      }
    }
    
    #Prevent vectors with less than 3 values
    if(length(unmethSampsD) < 3){
      unmethSampsD <- names(sort(meth)[1:3])
    }
    
    #create vector of RNA IDs linked to DNA IDs
    RNAsamp <- clinical[which(clinical$Tumor.Tissue.DNA.HAND.ID %in% colnames(Poodles)), ]
    DNAsamp <- RNAsamp$Tumor.Tissue.DNA.HAND.ID
    RNAsamp <- sapply(strsplit(RNAsamp$Tumor.Tissue.RNA.HAND.ID, split="\\_|-"), function(x){x[[1]]})
    names(RNAsamp) <- DNAsamp
    RNAsamp2 <- sort(RNAsamp[RNAsamp %in% names(Exp)])
    
    #subset meth vector
    meth <- meth[names(RNAsamp2)]
    
    #Only get Exp samples that also have meth data
    Exp <- Exp[RNAsamp2]
    
    #add avg promoter meth (for all samples) to avgMeth
    avgMeth[m] <- mean(meth)
    names(avgMeth)[m] <- a
    
    Kendall <- Kendall(Exp, meth)
    LM[m] <- Kendall$tau
    names(LM)[m] <- a
    LMP[m] <- Kendall$sl
    names(LMP)[m] <- a
    
    ###2nd Wilcoxon test with unmeth samples (exp vs. enhancer meth)
    
    #subset and sort unmethSamples
    #unmethSampsD2 = unmeth samps (DNA IDs) that are in RNAsamp (with rejected samples removed)
    unmethSampsD2 <- unmethSampsD[unmethSampsD %in% Allsamps]
    unmethSampsR <- RNAsamp2[unmethSampsD2]
    names(unmethSampsR) <- NULL
    
    #Exp vector for only unmeth samples
    newExp <- Exp[sort(unmethSampsR)]
    
    #Create new meth vector of average methylation of Enhancer 
    #(average of all noodles spanning enhancer) for each sample
    q <- 1
    newMeth <- numeric()
    
    for (z in unmethSampsD2){
      newMeth[q] <- mean(enhancerList[[i]][[7]][,z])
      names(newMeth)[q] <- RNAsamp2[z]
      q <- q+1
    }
    
    #sort newMeth
    newMeth <- newMeth[sort(names(newMeth))]
    
    Kendall2 <- Kendall(newExp, newMeth)
    LM2[m] <- Kendall2$tau
    names(LM2)[m] <- a
    LM2P[m] <- Kendall2$sl
    names(LM2P)[m] <- a
    
    m <- m + 1
    
    if (enhancerList[[i]][[20]] == "X"){
      break
    }
  } 
  
  ###### Create table for local meth info
  LMtable <- data.frame(AvgProMethylation = avgMeth, PromoterTau = LM, PromoterP = LMP, EnhancerTau = LM2, EnhancerP=LM2P)
  
  #add LMtable and wilcoxGenes to enhancerList
  enhancerList[[i]][[13]] <- LMtable

} 

save("enhancerList","sortedEnhancer", file = "LocalMethKendall2.Rda" )
