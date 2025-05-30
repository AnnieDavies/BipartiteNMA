#number of arms per study
num_arms <- function(data, trials, M){
  arms <- c()
  for(i in 1:M){
    dfi <- subset(data, trial==trials[i])
    arms <- c(arms, nrow(dfi))
  }
  arms
}

# Get mu from binary
get_lor_var <- function(data){
  data$LOR <- NA
  data$variance <- NA

  for(i in 1:nrow(data)){
    r<-data$r[i]+0.5
    n<-data$n[i]+0.5
    data$LOR[i] <- log(r/(n-r))
    data$variance[i] <- (r*(1-(r/n)))^(-1)
  }
  
  data$weight <- 1/data$variance
  data
}


#block diagonal variance matrix Sigma (within trial)
get_Sig <- function(data, trials, M, n){
  Sdiags <- list()
  rowcol_names <- c()
  for(i in 1:M){
    dfi <- subset(data, trial==trials[i])
    ti <- dfi$treatment
    Si <- matrix(0, nrow = n[i]-1, ncol = n[i]-1)
    #reference treatment is the first one
    refvar <- dfi$variance[1]
    for(j in 1:(n[i]-1)){
      for(k in 1:(n[i]-1)){
        if(j==k){
          Si[j,k] <- refvar + dfi$variance[j+1]
        }else{
          #off diag
          Si[j,k] <- refvar
        }
      }
    }
    rownames(Si) <- colnames(Si) <- paste(trials[i], ti[1], ti[2:n[i]])
    rowcol_names <- c(rowcol_names, rownames(Si))
    Sdiags[[i]] <- Si
  }
  S <- bdiag(Sdiags)
  rownames(S) <- colnames(S) <- rowcol_names
  S
}

#Between trial covariance matrix Omega
get_Omega <- function(data, trials, M, arms, tau){
  Odiags <- list()
  rowcol_names <- c()
  for(i in 1:M){
    dfi <- subset(data, trial==trials[i])
    ti <- dfi$treatment
    Oi <- matrix(0, nrow = arms[i]-1, ncol = arms[i]-1)
    for(j in 1:(arms[i]-1)){
      for(k in 1:(arms[i]-1)){
        if(j==k){
          Oi[j,k] <- tau^2
        }else{
          #off diag
          Oi[j,k] <- (tau^2)/2
        }
      }
    }
    rownames(Oi) <- colnames(Oi) <- paste(trials[i], ti[1], ti[2:arms[i]])
    rowcol_names <- c(rowcol_names, rownames(Oi))
    Odiags[[i]] <- Oi
  }
  O <- bdiag(Odiags)
  rownames(O) <- colnames(O) <- rowcol_names
  O
}


## define main design matrix
design_matrix <- function(data, trials, trts, M, N, arms){
  Xis <- list()
  for(i in 1:M){
    dfi <- subset(data, trial==trials[i])
    ti <- dfi$treatment
    
    Xi <- matrix(0, nrow = arms[i]-1, ncol = N-1)
    if(ti[1]==trts[1]){ #if trial reference = global reference
      for(j in 2:arms[i]){
        cind <- which(trts == ti[j])
        Xi[j-1,cind-1] <- 1
      }
    }else{ #otherwise we need -1 in trial reference column
      refind <- which(trts == ti[1]) #column of reference treatment
      for(j in 2:arms[i]){
        Xi[j-1, refind-1] <- -1
        
        cind <- which(trts == ti[j])
        Xi[j-1,cind-1] <- 1
      }
    }
    rownames(Xi)<- paste(trials[i], ti[1], ti[2:arms[i]])
    Xis[[i]] <- Xi
  }
  X <- do.call(rbind, Xis)
  colnames(X)<-paste(trts[1],trts[2:N])
  X
}


## contrast projection matrix
contrast_matrix <- function(data, trials, M, arms){
  Cdiags <- list()
  row_names <- c()
  col_names <- c()
  for(i in 1:M){
    dfi <- subset(data, trial==trials[i])
    ti <- dfi$treatment
    Ci <- get_CN(arms[i])
    rownames(Ci) <- paste(trials[i], ti[1], ti[2:arms[i]])
    colnames(Ci) <- paste(trials[i], ti[1:arms[i]])
    row_names <- c(row_names, rownames(Ci))
    col_names <- c(col_names, colnames(Ci))
    Cdiags[[i]]<-Ci
  }
  C <- as.matrix(bdiag(Cdiags))
  rownames(C) <- row_names
  colnames(C) <- col_names
  C
}

get_CN <- function(N){
  CN <- cbind(-rep(1,N-1), diag(N-1))
  CN
}


## define consistency eqs design matrix
consistency_matrix <- function(N, trts){
  D <- matrix(0, nrow=N*(N-1)/2, ncol = N-1)
  rownames(D)<-rep("name",nrow(D))
  colnames(D)<-rep("name",ncol(D))
  rind <- 0
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      rind <- rind + 1
      t1 <- trts[i]
      t2 <- trts[j]
      
      if(t1==trts[1]){
        cind <- which(trts==t2)
        D[rind, cind-1] <- 1
        colnames(D)[rind] <- paste(t1, t2)
      }else{
        refind <- which(trts==t1)
        D[rind, refind-1] <- -1
        cind <- which(trts==t2)
        D[rind, cind-1] <- 1
      }
      rownames(D)[rind] <- paste(t1, t2)
    }
  }
  D
}
