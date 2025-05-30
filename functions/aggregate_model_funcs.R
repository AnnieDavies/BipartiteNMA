#adjusted weights (reduce weights method)
get_Wadj <- function(data, trials, M, n, tau){
  rowcol_names <- c()
  Wdiags <- list()
  for(i in 1:M){
    dfi <- subset(data, trial==trials[i])
    ti <- dfi$treatment
    
    Vi <- var_eff_i(i, data, trials, n, tau) #effective variance matrix
    Lip <- L_i_plus(Vi) #Li+ in terms of Vi
    Li <- L_pinv(Lip) #take the pseudo-inverse
    Wi <- W_from_L(Li)
    Wdiags[[i]] <- Wi
    
    pairs <- c()
    for(j in 1:(length(ti)-1)){
      for(k in (j+1):length(ti)){
        pairs <- c(pairs, paste(ti[j], ti[k]))
      }
    }
    rownames(Wi) <- colnames(Wi) <- paste(trials[i], pairs)
    rowcol_names <- c(rowcol_names, rownames(Wi))
    
  }
  
  W <- bdiag(Wdiags)
  rownames(W) <- colnames(W) <- rowcol_names
  
  W
}

#"We add this estimate of τ2 to the observed sampling variance of each single 
# comparison in the network (before appropriately reducing the weights for multi-arm studies)"
var_eff_i <- function(i, data, trials, n, tau){

  dfi <- subset(data, trial==trials[i])
  Vi <- matrix(0, nrow = n[i], ncol = n[i])
  
  for(j in 1:n[i]){
    for(k in 1:n[i]){
      if(k==j){
        Vi[j,k] <- 0
      }else{
        Vi[j,k] <- dfi$variance[j]+dfi$variance[k]+tau^2
      }
    }
  }
  Vi
  
}

#Laplacian pseudo-inverse in terms of Vi
L_i_plus <- function(Vi){
  ni <- nrow(Vi)
  Oi <- matrix(1, nrow = ni, ncol = ni)
  term2 <- Vi%*%Oi + Oi%*%Vi
  term3 <- Oi%*%Vi%*%Oi
  Lip <- -0.5*(Vi - term2/ni + term3/(ni^2))
  Lip
}

#pseudo inverse of laplacian
L_pinv <- function(Li){
  ni <- nrow(Li)
  Oi <- matrix(1, nrow = ni, ncol = ni)
  pinv <- solve(Li - Oi/ni)+Oi/ni
  pinv
}

#Get weights from Laplacian matrix (-ve off diag entries)
W_from_L <- function(Li){
  mi <- nrow(Li)
  qi <- mi*(mi-1)/2
  
  Wdiags <- list()
  ind <- 1
  for(j in 1:(mi-1)){
    for(k in (j+1):mi){
      Wdiags[[ind]] <- -Li[j,k]
      ind <- ind+1
    }
  }
  Wi <- bdiag(Wdiags)
  Wi
}

#aggregate weight matrix
get_W_agg <- function(Wadj, trts, N, as_matrix=TRUE){
  rowcol_names <- c()
  Wdiag <- list()
  ind<-1
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      lab <- paste(trts[i], trts[j])
      t1t2 <- paste0(trts[i], " ", trts[j], "$")
      inds <-  grep(t1t2, colnames(Wadj))
      if(length(inds)>1){
        rowcol_names <- c(rowcol_names, lab)
        Wdiag[[ind]] <- sum(diag(Wadj)[inds])
        ind <- ind + 1
      }else if(length(inds)==1){
        rowcol_names <- c(rowcol_names, lab)
        Wdiag[[ind]] <- Wadj[inds,inds]
        ind <- ind + 1
      }
    }
  }
  Wa <- bdiag(Wdiag)
  rownames(Wa) <- colnames(Wa) <- rowcol_names
  Wa
}

#aggregate oriented edge-vertex incidence matrix
get_B_agg <- function(Wa, trts, N){
  K <- nrow(Wa)
  B <- matrix(0, nrow = K, ncol = N)
  for(i in 1:K){
    t1 <- strsplit(colnames(Wa)[i], " ")[[1]][1]
    t2 <- strsplit(colnames(Wa)[i], " ")[[1]][2]
    ind1 <- which(trts==t1)
    ind2 <- which(trts==t2)
    B[i, ind1] <- -1
    B[i, ind2] <- 1
  }
  rownames(B) <- colnames(Wa)
  colnames(B) <- trts
  B
}
