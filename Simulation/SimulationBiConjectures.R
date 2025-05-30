library(readxl)
library(igraph)
library(Matrix) #for bdiag
library(dplyr)
library(tidyverse)
library(pracma)
library(doParallel)
library(foreach)
library(parallel) 

## Load functions
path_func <-".../functions"
files <- list.files(path_func, full.names = TRUE)
lapply(files, function(f) {
  message("Sourcing: ", f)
  source(f)
})


## Set seed
set.seed(24601)

## Number of repetitions
REP <- 10000

## Dataframe to save results
cnames <- c("N", "M", "K_samp", "K_final", "s", "Max_diff_conj1", "Max_diff_conj2", 
            "Bi_density", "Uni_density", "Uni_radius", "Uni_distance", "Uni_edges",
            "Top_degree_mean", "Top_degree_min", "Top_degree_max", 
            "Bot_degree_mean", "Bot_degree_min", "Bot_degree_max",
            "Uni_degree_mean", "Uni_degree_min", "Uni_degree_max")
df <- as.data.frame(matrix(NA, nrow = REP, ncol = length(cnames) ))
colnames(df) <- cnames

## Simulation-------------------------------------------------------------------
for(i in 1:REP){
  print(i)
  
  #Sample NMA------------------
  N <- sample(3:50, size = 1) #sample number of treatments (bottom: type = TRUE)
  M <- sample(2:200, size =1) #sample number of trials (top: type = FALSE)
  sd <- (N*M-2*M)/2 
  K<-N*M+1 #number of arms
  while(K>N*M){
    K<-2*M+round(abs(rnorm(1, mean=0, sd = sd)))
  }
  s <- runif(1, min=0.5, max=2) #sample SD of variance distribution
  
  df$N[i] <- N
  df$M[i] <- M
  df$K_samp[i] <- K
  df$s[i] <- s
  
  g <- generate_NMA_bi(N, M, K, max_reps = 100, vardist = "halfnorm", varscale = s)
  narms <- ecount(g)
  df$K_final[i] <- narms
  
  # Define data---------------------
  data <- igraph::as_data_frame(g, what = "edges")
  colnames(data) <- c("trial", "treatment", "weight")
  data <- data[order(data$treatment), ]
  data <- data[order(data$trial), ]
  data$variance <- 1/data$weight
  
  # Extract variables
  # trials 
  trials <- unique(data$trial)
  # treatments
  trts <- unique(data$treatment)
  trts <- sort(trts) 
  
  
  ## Conjecture 1-----------------------
  # Arm level hat matrix
  arms <- num_arms(data, trials, M)
  Sigma <- get_Sig(data, trials, M, arms) #Within study variance matrix
  W <- solve(Sigma)
  X <- design_matrix(data, trials, trts, M, N, arms)
  C <- contrast_matrix(data, trials, M, arms)
  H_trial <- solve(t(X)%*%W%*%X)%*%t(X)%*%W
  H_arm <- H_trial%*%C
  
  # Matrix of edge currents
  Bbi <- get_bi_incidence(data, trials, trts, oriented = TRUE)
  R <- get_R(data)
  Rinv <- matrix(0, nrow=nrow(R), ncol=ncol(R))
  diag(Rinv)<-1/diag(R)
  rownames(Rinv) <- colnames(Rinv) <- rownames(R)
  Jt <- get_Jt(trts, trials)
  L <- t(Bbi)%*%Rinv%*%Bbi
  It <- Jt %*% pinv(L) %*% t(Bbi) %*% Rinv
  
  # Max difference
  df$Max_diff_conj1[i] <- max(abs(H_arm-It))
  
  ## Conjecture 2------------------------
  # T (transition matrix on unipartite graph)
  Wadj <- get_Wadj(data, trials, M, arms, tau=0)
  Wa <- get_W_agg(Wadj, trts, N)
  T.agg <- agg_trans_matrix(Wa, N, trts)
  
  # P tilde (renormalised 2-step transition matrix)
  B <- bi_adj_matrix(data, M, N, trials, trts, weighted = TRUE)
  P.up <- t(B)/rowSums(t(B)) #upstream: from treatments to trials
  P.down <- B/rowSums(B)#downstream: from trials to treatments
  P <- P.up%*%P.down 
  P.rn <- renorm_P(P)
  
  # Max difference
  df$Max_diff_conj2[i] <- max(abs(T.agg-P.rn))
  
  ## Network characteristics------------------------
  df$Bi_density[i] <- narms/(N*M)
  
  #Top nodes
  top_nodes <- V(g)[type==FALSE]
  top_degrees <- degree(g, top_nodes)
  df$Top_degree_mean[i] <- mean(top_degrees)
  df$Top_degree_min[i] <- min(top_degrees)
  df$Top_degree_max[i] <- max(top_degrees)
  
  #Bottom degrees
  bot_nodes <- V(g)[type==TRUE]
  bot_degrees <- degree(g, bot_nodes)
  df$Bot_degree_mean[i] <- mean(bot_degrees)
  df$Bot_degree_min[i] <- min(bot_degrees)
  df$Bot_degree_max[i] <- max(bot_degrees)
  
  # bottom projection
  g.bot <- bipartite_projection(g, multiplicity = FALSE, which = "true")
  
  df$Uni_density[i] <- edge_density(g.bot)
  df$Uni_radius[i] <- radius(g.bot)
  df$Uni_distance[i] <- mean_distance(g.bot)
  
  #Unipartite degree
  uni_degrees <- degree(g.bot, V(g.bot))
  df$Uni_degree_mean[i] <- mean(uni_degrees)
  df$Uni_degree_min[i] <- min(uni_degrees)
  df$Uni_degree_max[i] <- max(uni_degrees)
  
  ##number of edges (unipartite)
  df$Uni_edges[i] <- ecount(g.bot)
}

# Conjectures: maximum differences
max(df$Max_diff_conj1, na.rm=TRUE) #6.027618e-11
max(df$Max_diff_conj2, na.rm=TRUE) #1.147971e-13

# Save data
rnames <- c("Minimum", "Maximum", "Mean", "SD")
df.res <- data.frame(
  Minimum = apply(df, 2, min, na.rm = TRUE),
  Maximum = apply(df, 2, max, na.rm = TRUE),
  Mean = apply(df, 2, mean, na.rm = TRUE),
  SD = apply(df, 2, sd, na.rm = TRUE)
)

df.res <- t(df.res)  # Transpose to get the desired format
df.res <- as.data.frame(df.res)
rownames(df.res)<-rnames
