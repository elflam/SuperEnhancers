###3/24/2017
###Emily Flam
###Same as VisualsAltered but modified to only include T and N samples

###Functions
#serially merge multiple matrices into single matrix 
expmat <- function(x, y){
  df <- merge(x, y, by= "row.names", all.x= F, all.y= F, sort=F)
  rownames(df) <- df$Row.names
  df$Row.names <- NULL
  return(df)
}

###Main Code
setwd("~/Desktop/Super-Enhancers/Code Skeleton")
load("LocalMethKendall2.Rda")

#Load noodle report frame
#huge.report.frame contains coordinate and meth ratio
load("/Volumes/Seagate Backup Plus Drive/Noodles.huge.report.frame.Rda")

library(gplots)

#Create Folder in Super-Enhancers to hold all files from run; set as working directory
dir.create(paste0("~/Desktop/Super-Enhancers/",Sys.Date(),"Run",sep=""))

for (i in 1:length(enhancerList)){
  k <- i
  
  setwd(paste0("~/Desktop/Super-Enhancers/",Sys.Date(),"Run"))
  
  pdf(file=paste0(enhancerList[[i]][[3]],"_", enhancerList[[i]][[4]], "-", enhancerList[[i]][[5]]))
  
  #create label for pdf 
  label <- paste(enhancerList[[i]][[3]], ":", enhancerList[[i]][[4]], "-", enhancerList[[i]][[5]])
  textplot(label)
  
  #Collect noodles (huge.report.frame and noodles.C.methylation) that overlap with enhancer
  noodleOverlap <- huge.report.frame[which(huge.report.frame$end >= enhancerList[[i]][[4]] & huge.report.frame$start <= enhancerList[[i]][[5]] & huge.report.frame$chr == enhancerList[[i]][[3]]), ]
  
  #create plot with both normal & tumor meth ratios
  Coordinate <- cbind(noodleOverlap$start, noodleOverlap$start)
  Methylation <- cbind(noodleOverlap$nor.ratio, noodleOverlap$tmr.ratio)
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
  matplot(Coordinate,Methylation, type="l", ylim=c(0,1))
  points(c(noodleOverlap[enhancerList[[i]][[6]], ]$start), rep(1, length(enhancerList[[i]][[6]])), col = "blue", pch = 8)
  legend('topright', inset=c(-0.25,0), pch='-',col=c('red','black'),legend=c('T','N'))
  
  #modify for 1-row matrices
  if (enhancerList[[k]][[20]] == "X"){
    
    heatmap.2(enhancerList[[i]][[1]], Colv=F, Rowv=F, trace="none", scale='row', col=redgreen, ColSideColors = enhancerList[[i]][[10]], hclust=function(x) hclust(x,method="complete"), distfun=function(x) as.dist((1-cor(t(x)))/2))
    
    #number of waterfall plots to make
    wf <- 1
    
  } else{
    #create heatmap
    heatmap.2(enhancerList[[i]][[1]], Colv=F, Rowv=F, trace="none", scale='row', col=redgreen, ColSideColors = enhancerList[[i]][[10]], hclust=function(x) hclust(x,method="complete"), distfun=function(x) as.dist((1-cor(t(x)))/2))
    
    #number of waterfall plots to make
    if (nrow(enhancerList[[k]][[12]]) >= 5){
      wf <- 5
    }else{
      wf <- nrow(enhancerList[[k]][[12]])
    }
  }
  
  #create waterfall plots
  for (i in 1:wf){
    samples <- names(enhancerList[[k]][[2]])[order(enhancerList[[k]][[2]], enhancerList[[k]][[1]][rownames(enhancerList[[k]][[12]][i,]),])]
    x = as.numeric(enhancerList[[k]][[1]][rownames(enhancerList[[k]][[12]][i,]), samples])
    names(x) <- samples
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
    barplot(x, las=2, col = enhancerList[[k]][[10]][samples])
    legend(x="topright", inset=c(-0.25,0), legend = unique(enhancerList[[k]][[2]][order(enhancerList[[k]][[2]])]),
    fill = unique(enhancerList[[k]][[10]][samples]), title = rownames(enhancerList[[k]][[12]][i,]))
    }
  
  #Combine DEStats and Wilcoxon table into one matrix; add extra column with enhancer info
  enhancerList[[k]][[14]] <- expmat(enhancerList[[k]][[12]],enhancerList[[k]][[13]])
  enhancerList[[k]][[14]][,13] <- c(rep(paste0(enhancerList[[k]][[3]], ":", enhancerList[[k]][[4]], "-", enhancerList[[k]][[5]]),nrow(enhancerList[[k]][[14]])))
  colnames(enhancerList[[k]][[14]])[13] <- "Enhancer"
  enhancerList[[k]][[14]][,14] <- rownames(enhancerList[[k]][[14]])
  colnames(enhancerList[[k]][[14]])[14] <- "Gene"
  rownames(enhancerList[[k]][[14]]) <- NULL
    
  #print gene table
  textplot(enhancerList[[k]][[12]])
  textplot(enhancerList[[k]][[13]])
  
  dev.off()
  
  #Add Genes to big gene table; use first enhancer to set colnames and fill thereafter
  if(k==1){
    BigGeneTable <- enhancerList[[k]][[14]]
  }else{
    BigGeneTable <- rbind(BigGeneTable, enhancerList[[k]][[14]])
  }
}  

write.table(BigGeneTable, file=paste0("~/Desktop/Super-Enhancers/",Sys.Date(),"Run/BigGeneTable.txt"), sep="\t")
write.table(sortedEnhancer, file=paste0("~/Desktop/Super-Enhancers/",Sys.Date(),"Run/sortedEnhancersTvN.txt"),sep = ",")
  