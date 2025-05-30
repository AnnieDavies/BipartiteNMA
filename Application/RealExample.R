## Annabel L Davies 2025

library(dbplyr)
library(multinma)
library(tidyverse)
library(igraph)
library(netmeta)
library(Matrix) 
library(pracma)
library(xtable)
library(ggraph)
library(ggplot2)
options(mc.cores = parallel::detectCores())

# Load functions
path_func <- "...\\functions"
files <- list.files(path_func, full.names = TRUE)
lapply(files, function(f) {
  message("Sourcing: ", f)
  source(f)
})

# Define data-------------------------------------------------------------------
# Read in plaque psoriasis data from multinma package
data <- plaque_psoriasis_agd #aggregate data
data_ipd <- plaque_psoriasis_ipd #ipd

#Convert IPD trials to aggregate data
ipd_arms <- group_by(plaque_psoriasis_ipd, studyc, trtc) %>% 
  summarise(pasi75_r = sum(pasi75),
            pasi75_n = n())

data <- bind_rows(
  select(data, !! colnames(ipd_arms)),
  ipd_arms
)

colnames(data) <- c("trial", "treatment", "r", "n")
data <- data[order(data$treatment), ]
data <- data[order(data$trial), ]


trials <- unique(data$trial)
trts <- unique(data$treatment)
trts <- sort(trts)


M <- length(trials) #Number of trials
N <- length(trts) #Number of treatments


# Calculate arm level means and variances
data <- get_lor_var(data)


# Arm and trial level matrices--------------------------------------------------
mu <- data$LOR #arm level means
arms <- num_arms(data, trials, M)
C <- contrast_matrix(data, trials, M, arms)
y <- C%*%mu #constrast level data

Sigma <- get_Sig(data, trials, M, arms) #Within study variance matrix
W <- solve(Sigma)
X <- design_matrix(data, trials, trts, M, N, arms) #Design matrix nrow = sum_i arms[i]-1, ncol = N-1

#Trial-level hat matrix
H_trial <- solve(t(X)%*%W%*%X)%*%t(X)%*%W

#Arm-level hat matrix
H_arm <- H_trial%*%C

# Graph matrices----------------------------------------------------------------
#biadjacency matrix
B <- bi_adj_matrix(data, M, N, trials, trts, weighted = TRUE)
#bipartite incidence matrix
Bbi <- get_bi_incidence(data, trials, trts, oriented = TRUE)

# Electric network matrices-----------------------------------------------------
R <- get_R(data) #resistances
Rinv <- solve(R)
Jt <- get_Jt(trts, trials) #nodal currents
L <- t(Bbi)%*%Rinv%*%Bbi

It <- Jt %*% pinv(L) %*% t(Bbi) %*% Rinv #edge currents


########### Conjecture 1 ################
max_diff1 <- max(abs(H_arm-It))


# Aggregate matrices------------------------------------------------------------
# Adjusted weight matrix
Wadj <- get_Wadj(data, trials, M, arms, tau=0)
# Aggregate weight matrix
Wa <- get_W_agg(Wadj, trts, N)

# Incidence matrix of unipartite graph
Ba <- get_B_agg(Wa, trts, N)

La <- t(Ba)%*%Wa%*%Ba
Lpa <- L_pinv(La)
CN <- get_CN(N)
rownames(CN) <- paste(trts[1], trts[2:N])
colnames(CN) <- trts

## Aggregate Hat matrix
H_agg <- CN%*%Lpa%*%t(Ba)%*%Wa

# Transition matrices-----------------------------------------------------------
## Bipartite RW
# each element ij is the probability of going from i->j
P.up <- t(B)/rowSums(t(B)) #upstream: from treatments to trials
P.down <- B/rowSums(B)#downstream: from trials to treatments

## Two-step transition matrix
P <- P.up%*%P.down 

## Aggregate transition matrix
T.agg <- agg_trans_matrix(Wa, N, trts)

## Renormalise P
P.rn <- renorm_P(P) ##we need P.rn = T.agg


############ Conjecture 2 ############
max_diff2 <- max(abs(T.agg-P.rn))



# Plots: evidence flow----------------------------------------------------------

#Unipartite flow----------------------------------------------------------------

## Consistency equations matrix (to obtain all comparisons)
D <- consistency_matrix(N, trts)

H_agg_full<-D%*%H_agg

F26.uni <- get_uniFlow(D%*%H_agg, source = trts[2], sink = trts[6], N)
g.uniflow26 <- graph_from_adjacency_matrix(F26.uni, mode = "directed", weighted = TRUE, diag = FALSE)

uf <- plot_uniflow_gg(g.uniflow26, source = trts[2], sink = trts[6], layout = "circle",
                      vsize=12, vtxt=8, etxt = 6, ashift=7, dodge = 3, push=7, 
                      border_margin = unit(c(1, 1, 1, 1), "cm"), col.edge = "#707070", 
                      emin=0.2, emax=2.5)


#Bipartite flow-----------------------------------------------------------------
H_arm_full<-D%*%H_arm

F26.bi <- get_biFlow(D%*%H_arm, source = trts[2], sink = trts[6], N, M, trials, trts, order = 2) #treatment then trial

# Create weighted directed bipartite graph
g.biflow26 <- graph_from_adjacency_matrix(F26.bi, mode = "directed", weighted = TRUE, diag = FALSE)
V(g.biflow26)$type <- c(rep(TRUE, N), rep(FALSE, M)) #treatment then trial

bf <- plot_biflow_gg(g.biflow26, trts[2], trts[6], layout = "bipartite", vsize=12, vtxt=8, ashift=7, etxt = 6, col.edge = "#707070",
                     dodge = 3, push = -2, emin=0.2, emax=2.5, border_margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

