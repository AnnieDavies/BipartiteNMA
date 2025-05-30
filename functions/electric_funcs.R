get_R <- function(data){
  R <- diag(1/data$weight)
  colnames(R) <- rownames(R) <- paste(data$trial, data$treatment)
  R
}

#edge-vertex incidence of bipartite graph
get_bi_incidence <- function(data, trials, trts, oriented = TRUE){
  #N=number of treatments, M=number of trials, K=number of arms (edges)
  K<-nrow(data)
  N<-length(trts)
  M<-length(trials)
  
  Bbi <- matrix(0, nrow=K, ncol = N+M)
  rownames(Bbi)<-paste(data$trial, data$treatment)
  colnames(Bbi)<-c(trials, trts)
  
  for(i in 1:K){
    source_ind <- which(colnames(Bbi)==data$treatment[i])
    sink_ind <- which(colnames(Bbi)==data$trial[i])
    Bbi[i, source_ind] <- 1
    if(oriented){
      Bbi[i, sink_ind] <- -1
    }else{
      Bbi[i, sink_ind] <- 1
    }
  }
  Bbi 
}


get_Jt <- function(trts, trials){
  N <- length(trts)
  M <- length(trials)
  
  CN <- get_CN(N)
  Jt <- cbind(matrix(0, nrow = N-1, ncol = M), CN)
  rownames(Jt) <- paste(trts[1], trts[2:N])
  colnames(Jt)<-c(trials, trts)
  Jt
}