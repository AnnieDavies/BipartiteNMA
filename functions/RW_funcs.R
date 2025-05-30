full_bi_trans_matrix <- function(P.up, P.down, N, M, nodes){
  # Create zero matrices for the block diagonals
  zero_top <- matrix(0, nrow = N, ncol = N)  # Zero matrix for top-left block
  zero_bot <- matrix(0, nrow = M, ncol = M)  # Zero matrix for bottom-right block
  
  # Combine the matrices to form the desired block matrix
  # The order is:
  # [zero_top | U]
  # [D        | zero_bot]
  Trans <- rbind(
    cbind(zero_top, P.up), 
    cbind(P.down, zero_bot))
  rownames(Trans)<-colnames(Trans)<-nodes
  
  Trans
}

agg_trans_matrix <- function(Wa, N, trts){
  T.agg <- matrix(0, nrow = N, ncol = N)
  rownames(T.agg) <- colnames(T.agg) <- trts
  
  W.vec <- diag(Wa)
  
  for(j in 1:N){
    
    inds <- grep(paste0(trts[j], "( |$)"), names(W.vec))
    names(W.vec[inds])
    if(length(inds)>1){
      denom <- sum(W.vec[inds])
    }else if(length(inds)==1){
      denom <- W.vec[inds]
    }else{
      denom <- 0
    }
    
    for(k in 1:N){
      if(j==k){
        T.agg[j,k] <- 0
      }else{
        pattern <- paste0("(", trts[j], " ", trts[k], "$|", trts[k], " ", trts[j], "$)")
        inds2 <- grep(pattern, names(W.vec))
        if(length(inds2)>0){
          T.agg[j,k] <- W.vec[inds2]/denom
        }else{
          T.agg[j,k] <- 0
        }
      }
    }
  }
  T.agg
}

renorm_P <- function(P){
  diag(P)<-0
  P <- t(apply(P, 1, function(row) row / sum(row)))
  P
}

sim_RW_undir <- function(Trans, nodes, source_node, sink_node, n_sim, progress){
  # Initialize the count matrix to store net crossings
  tot_net_cross <- matrix(0, nrow=length(nodes), ncol=length(nodes))
  rownames(tot_net_cross) <- nodes
  colnames(tot_net_cross) <- nodes
  
  for(sim in 1:n_sim){
    if(sim %% progress==0){
      print(sim)
    }
    
    # Initialize the count matrix 
    net_cross <- matrix(0, nrow=length(nodes), ncol=length(nodes))
    
    # Start from the source_node node
    current_node <- source_node
    
    while (current_node != sink_node) {
      # Find index of current node 
      current_index <- which(nodes == current_node)
      
      # Choose the next node based on the transition probabilities
      next_node <- sample(nodes, 1, prob=Trans[current_index, ])
      
      # Determine the next node index
      next_index <- which(nodes == next_node)
      
      # Update net crossings for this simulation
      net_cross[current_index, next_index] <- net_cross[current_index, next_index] + 1
      net_cross[next_index, current_index] <- net_cross[next_index, current_index] - 1
      
      # Move to the next node
      current_node <- next_node
    }
    
    # Accumulate the results
    tot_net_cross <- tot_net_cross + net_cross
  }
  
  # Calculate the average crossings
  av_net_cross <- tot_net_cross / n_sim
  av_net_cross[av_net_cross < 0] <- 0
  
  av_net_cross
}

