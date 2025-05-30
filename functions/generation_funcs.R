generate_NMA_bi <- function(N, M, K, max_reps = 50, vardist = "halfnorm", varscale = 2){
  
  g <- generate_graph(N, M, K, max_reps)
  
  narms <- ecount(g)
  if(vardist=="halfnorm"){
    E(g)$weight <- 1/abs(rnorm(n = narms, mean = 0, sd = varscale))
  }else if(vardist=="unif"){
    E(g)$weight <- 1/runif(narms, min = varscale[1], max = varscale[2])
  }else{
    print("Error: vardist must be 'halfnorm' or 'unif'")
  }
  
  g
}

generate_graph <- function(N, M, K, max_reps){
  g <- sample_bipartite(n1 = M, n2 = N, type = "gnm", m = K, directed = FALSE)
  
  #assign names to nodes
  V(g)$name[V(g)$type] <- paste0("T", 1:N)  # For bottom nodes (type = TRUE)
  V(g)$name[!V(g)$type] <- paste0("S", 1:M) # For top nodes (type = FALSE)
  
  # Ensure that every trial contains >1 arm-------------------------------------
  while (any(degree(g, !V(g)$type) < 2)) {
    # Find top nodes with degree < 2
    top_nodes <- V(g)[!V(g)$type]
    edit.trials <- top_nodes[degree(g, top_nodes) < 2] 
    
    # Step 3: Add edges for underconnected top nodes
    for (v in edit.trials) {
      # Get potential treatments that are not yet in this trial
      bottom_nodes <- V(g)[V(g)$type]  # Bottom nodes (type = TRUE)
      not_adjacent <- sapply(bottom_nodes, function(x) !are_adjacent(g, v, x)) # Nodes not adjacent to v
      poss.trts <- bottom_nodes[not_adjacent] 
      
      if (length(poss.trts) >= 2 - degree(g, v)) {
        # Randomly select bottom nodes to form edges
        new_edges <- sample(poss.trts, size = 2 - degree(g, v), replace = FALSE)
        
        # Add the new edges
        g <- add_edges(g, c(rbind(v, new_edges)))
      }
    }
  }
  
  #Ensure the treatments are connected------------------------------------------
  
  #count_components 
  n.comp <- count_components(g)
  
  rep<-1
  while(n.comp > 1 && rep<max_reps){
    components <- components(g)
    comp_membership <- components$membership
    
    #list of the different components
    unique_comps <- unique(comp_membership)
    
    #randomly sample 2 components
    rand_comps <- sample(unique_comps, 2, replace = FALSE)
    
    # Sample a bottom node in component 1 (if it has more than one node)
    comp1 <- V(g)[which(comp_membership == rand_comps[1])]
    comp1 <- comp1[type==TRUE] #bottom nodes are type TRUE
    if(length(comp1) > 1){
      node1 <- sample(comp1, size = 1)
    }else{
      node1 <- comp1
    }
    
    # Sample a bottom node in component 2 if it has more than one node
    comp2 <- V(g)[which(comp_membership == rand_comps[2])]
    comp2 <- comp2[type==TRUE] #bottom nodes are type TRUE
    if(length(comp2) > 1){
      node2 <- sample(comp2, size = 1)
    }else{
      node2 <- comp2
    }
    
    #All trials
    top_nodes <- V(g)[!V(g)$type]
    # Find trials connected to node1 
    top1 <- top_nodes[sapply(top_nodes, function(node) {are_adjacent(g, node, node1)})]
    # Find trials connected to node2
    top2 <- top_nodes[sapply(top_nodes, function(node) {are_adjacent(g, node, node2)})]
    tops <- c(top1,top2)
    
    if(length(tops)>0){
      
      #sample a trial
      if(length(tops)>1){
        trial <- sample(tops, size = 1)
      }else{
        trial <- tops
      }
      
      #find out if it is connected to 1 or 2 and then add an edge to the other
      if(are_adjacent(g, trial, node1)){
        #trial already contains node1 so add an arm to node2
        g <- add_edges(g, c(trial, node2))
      }else if(are_adjacent(g, trial, node2)){
        #trial already contains node2 so add an arm to node1
        g <- add_edges(g, c(trial, node1))
      }else{
        print("Error: sampled trial is not connected to either node")
      }
    }else{
      #if there are no trials connected to either then add two arms (edges) to a trial (selected at random)
      trial <- sample(top_nodes, size = 1)
      g <- add_edges(g, c(trial, node1, trial, node2))
    }
    
    #count_components (we want the bottom projection to have one connected component)
    n.comp <- count_components(g)
    rep<-rep+1
  }
  
  g
  
}