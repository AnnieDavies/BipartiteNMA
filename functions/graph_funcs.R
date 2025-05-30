bi_adj_matrix <- function(data, M, N, trials, trts, weighted = TRUE){
  B <- matrix(0, nrow = M, ncol = N)
  rownames(B) <- trials
  colnames(B) <- trts
  for(i in 1:M){
    dfi <- subset(data, trial==trials[i])
    ti <- dfi$treatment
    wi <- dfi$weight
    
    for(j in 1:length(ti)){
      rind <- which(trts==ti[j])
      if(weighted){
        B[i, rind] <- wi[j] #weighted adjacency matrix
      }else{
        B[i, rind] <- 1
      }
    }
  } 
  B
}

get_biFlow <- function(H_bi, source, sink, N, M, trials, trts, order){
  
  Frow <- which(rownames(H_bi)==paste(source, sink))
  Fmat <- matrix(0, nrow = N+M, ncol = N+M)
  
  # Assign row and column names
  if(order == 1){#trial then treatment
    rownames(Fmat) <- c(trials, trts)
    colnames(Fmat) <- c(trials, trts)
  }else if(order==2){ #treatment then trial
    rownames(Fmat) <- c(trts, trials)
    colnames(Fmat) <- c(trts, trials)
  }
  
  # Populate the adjacency matrix based on the Frow row of the H matrix
  for (i in 1:ncol(H_bi)) {
    label <- unlist(strsplit(colnames(H_bi)[i], " "))
    trial_idx <- which(rownames(Fmat) == label[1])
    trt_idx <- which(colnames(Fmat) == label[2])
    
    if(H_bi[Frow, i] > 0){
      Fmat[trial_idx, trt_idx] <- H_bi[Frow, i]
    } else {
      Fmat[trt_idx, trial_idx] <- -H_bi[Frow, i]
    }
  }
  Fmat
}


#Assign aggregate edge weights to the unipartite graph
assign_agg_weights <- function(g, Wa){
  edges <- as_edgelist(g)
  edge_labels <- apply(edges, 1, function(x) paste(sort(x), collapse = " "))
  edge_weights <- sapply(edge_labels, function(e) {
    if (e %in% rownames(Wa)) {
      Wa[e, e]
    } else {
      NA  
    }
  })
  E(g)$weight <- edge_weights
  g
}

get_uniFlow <- function(Ha, source, sink, N){
  
  Frow <- which(rownames(Ha)==paste(source, sink))
  Fmat <- matrix(0, nrow = N, ncol = N)
  
  # Assign row and column names
  rownames(Fmat) <- colnames(Fmat) <- trts
  
  # Populate the adjacency matrix based on the Frow row of the H matrix
  for (i in 1:ncol(Ha)) {
    label <- unlist(strsplit(colnames(Ha)[i], " "))
    t1 <- which(rownames(Fmat) == label[1])
    t2 <- which(colnames(Fmat) == label[2])
    
    if(Ha[Frow, i] > 0){
      Fmat[t1, t2] <- Ha[Frow, i]
    } else {
      Fmat[t2, t1] <- -Ha[Frow, i]
    }
  }
  Fmat
}


plot_biflow_gg <- function(g.flow, source, sink, layout = "sugiyama",
                           vsize=10, vtxt=6, etxt = 5, emin=0.2, emax=2,
                           ashift=5, dodge=3, push = 0, digits=3,
                           col.trt = "#F5CAC3", col.trial = "#84A59D", col.sourcesink = "#F28482",
                           col.vtxt = "black", col.etxt = "black", col.edge = "black",
                           efont = "sourcesanspro", vfont = "sourcesanspro",
                           border_margin = unit(c(0, 0, 0, 0), "cm")){
                           
  vertex_colors <- ifelse(V(g.flow)$type, col.trt, col.trial)
  vertex_colors[V(g.flow)$name %in% c(source, sink)] <- col.sourcesink  # color source and sink
  
  ggraph(g.flow, layout=layout) +
    geom_edge_link(aes(width = weight, label = round(weight,digits)),
                   arrow = arrow(),
                   angle_calc = 'along',
                   label_size = etxt,
                   label_colour = col.etxt,
                   color = col.edge,
                   family = efont,
                   label_dodge = unit(dodge, 'mm'),
                   label_push = unit(push, 'mm'),
                   end_cap = circle(ashift, "mm")) +
    geom_node_point(aes(color = I(vertex_colors)), size = vsize,
                    shape = I(ifelse(V(g.flow)$type, "circle", "square")))+
    geom_node_text(aes(label = name), color = col.vtxt, size = vtxt, family = vfont)+
    scale_edge_width(range = c(emin, emax)) +
    theme_void() +
    theme(legend.position = "none", plot.margin = border_margin)+
    coord_cartesian(clip = "off", expand = TRUE)
  
}


plot_uniflow_gg <- function(g.flow, source, sink, layout = "circle",
                           vsize=10, vtxt=6, etxt = 5, emin=0.2, emax=2,
                           ashift=5, dodge = 3, push=6, digits=3,
                           col.trt = "#F5CAC3", col.sourcesink = "#F28482",
                           col.vtxt = "black", col.etxt = "black", col.edge = "black",
                           efont = "sourcesanspro", vfont = "sourcesanspro",
                           border_margin = unit(c(0, 0, 0, 0), "cm")){

  vertex_colors <- rep(col.trt, length(V(g.flow)))
  vertex_colors[V(g.flow)$name %in% c(source, sink)] <- col.sourcesink #colour source and sink

  ggraph(g.flow, layout=layout) +
    geom_edge_link(aes(width = weight, label = round(weight,digits)),
                   arrow = arrow(),
                   angle_calc = 'along',
                   label_size = etxt,
                   label_colour = col.etxt,
                   color = col.edge,
                   family = efont,
                   label_dodge = unit(dodge, 'mm'),
                   label_push = unit(push, 'mm'),
                   end_cap = circle(ashift, "mm")) +
    geom_node_point(aes(color = I(vertex_colors)), size = vsize,
                    shape = "circle")+
    geom_node_text(aes(label = name), color = col.vtxt, size = vtxt, family = vfont)+
    scale_edge_width(range = c(emin, emax)) +
    theme_void() +
    theme(legend.position = "none", plot.margin = border_margin)+
    coord_cartesian(clip = "off", expand = TRUE)
}


plot_bi_gg <- function(g.bi, layout = "bipartite",
                           vsize=10, vtxt=6, etxt = 5, emin=0.2, emax=2,
                           endcap=5, startcap=5, digits=3, dodge = 3, push = 12,
                           col.trt = "#F5CAC3", col.trial = "#84A59D",
                           col.vtxt = "black", col.etxt = "black", col.edge = "black",
                           efont = "sourcesanspro", vfont = "sourcesanspro",
                          border_margin = unit(c(0, 0, 0, 0), "cm")){
  
  vertex_colors <- ifelse(V(g.bi)$type, col.trt, col.trial)
  ggraph(g.bi, layout=layout) +
    geom_edge_link(aes(width = weight, label = round(weight,digits)),
                   angle_calc = 'along',
                   label_size = etxt,
                   label_colour = col.etxt,
                   color = col.edge,
                   family = efont,
                   label_dodge = unit(dodge, 'mm'),
                   label_push = unit(push, 'mm'),
                   end_cap = circle(endcap, "mm"),
                   start_cap = circle(startcap, "mm")) +
    geom_node_point(aes(color = I(vertex_colors)), size = vsize,
                    shape = I(ifelse(V(g.bi)$type, "circle", "square")))+
    geom_node_text(aes(label = name), color = col.vtxt, size = vtxt, family = vfont)+
    scale_edge_width(range = c(emin, emax)) +
    theme_void() +
    theme(legend.position = "none", plot.margin = border_margin)+
    coord_cartesian(clip = "off", expand = TRUE)
}

plot_uni_gg <- function(g, layout = "bipartite",
                       vsize=10, vtxt=6, etxt = 5, emin=0.2, emax=2,
                       endcap=5, startcap=5, digits=3, dodge = 3, push = 12,
                       col.trt = "#F5CAC3", col.vtxt = "black", col.etxt = "black", col.edge = "black",
                       efont = "sourcesanspro", vfont = "sourcesanspro",
                       border_margin = unit(c(0, 0, 0, 0), "cm")){
  
  vertex_colors <- rep(col.trt, length(V(g)))
  ggraph(g, layout=layout) +
    geom_edge_link(aes(width = weight, label = round(weight,digits)),
                   angle_calc = 'along',
                   label_size = etxt,
                   label_colour = col.etxt,
                   color = col.edge,
                   family = efont,
                   label_dodge = unit(dodge, 'mm'),
                   label_push = unit(push, 'mm'),
                   end_cap = circle(endcap, "mm"),
                   start_cap = circle(startcap, "mm")) +
    geom_node_point(aes(color = I(vertex_colors)), size = vsize,
                    shape = "circle")+
    geom_node_text(aes(label = name), color = col.vtxt, size = vtxt, family = vfont)+
    scale_edge_width(range = c(emin, emax)) +
    theme_void() +
    theme(legend.position = "none", plot.margin = border_margin)+
    coord_cartesian(clip = "off", expand = TRUE)
}