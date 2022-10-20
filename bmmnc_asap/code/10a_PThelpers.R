library(Matrix)
library(SummarizedExperiment)
library(tidyverse)
library(umap)
library(edgeR)
library(FNN)
library(matrixStats)
library(igraph)
library(BuenColors)
library(ggrastr)
set.seed(1)
"%ni%" <- Negate("%in%")

getQuantiles <- function(x){
  trunc(rank(x))/length(x)
}

alignTrajectory <- function(df, trajectory, filter = 0.05, dof = 250, spar = 1){
  findClosest <- function(x, y, fitx, fity){
    distxy <- sqrt(rowSums(cbind((fitx - x)^2 + (fity - y)^2)))
    idxMin <- which.min(distxy)
    if(idxMin==1){
      idxMin <- idxMin + 1
    }else if(idxMin==length(fitx)){
      idxMin <- idxMin - 1
    }
    if(distxy[idxMin + 1]  < distxy[idxMin - 1]){
      diff <- 1
    }else{
      diff <- -1
    }
    data.frame(idx = idxMin, dist = distxy[idxMin], diff = diff)
  }
  dfAll <- data.frame()
  for(x in seq_along(trajectory)){
    #Subset
    dfx <- df[df$Group==trajectory[x],]
    #Mean Diff Filter
    xmean <- colMeans(dfx[,c(1,2)])
    diffx <- sqrt(colSums((t(dfx[,1:2]) - xmean)^2))
    dfx <- dfx[which(diffx <= quantile(diffx,1 - filter)),]
    #Get diff
    if(x!=length(trajectory)){
      xmean1 <- colMeans(df[df$Group==trajectory[x+1],c(1,2)])
      diffx1 <- sqrt(colSums((t(dfx[,1:2]) - xmean1)^2))
      dfx$time <- (1 - getQuantiles(diffx1)) + x
    }else{
      xmean1 <- colMeans(df[df$Group==trajectory[x-1],c(1,2)])
      diffx1 <- sqrt(colSums((t(dfx[,1:2]) - xmean1)^2))
      dfx$time <- getQuantiles(diffx1) + x
    }
    dfAll <- rbind(dfAll , dfx)
  }
  sx <- smooth.spline(dfAll$time, dfAll$x, df = dof, spar = spar)
  sy <- smooth.spline(dfAll$time, dfAll$y, df = dof, spar = spar)
  dfFit <- data.frame(x = sx[[2]], y = sy[[2]], t = seq_along(sy[[2]]))
  dfTrajectory <- df[df$Group %in% trajectory,]
  dfTime <- lapply(seq_len(nrow(dfTrajectory)), function(x){
    findClosest(dfTrajectory[x,1],dfTrajectory[x,2], dfFit[,1],dfFit[,2])
  }) %>% Reduce("rbind",.)
  dfTime$distQ <- getQuantiles(dfTime$dist)
  dfTrajectory$pseudotime <- 100*getQuantiles(dfTime$idx + matrixStats::rowProds(as.matrix(dfTime[,c("diff","distQ")])))
  
  out <- list(trajectory=dfTrajectory, fitTrajectory = dfFit)
}
