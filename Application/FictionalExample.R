## Annabel L Davies 2025

## Load libraries
library(readxl)
library(igraph)
library(ggraph)
library(ggplot2)
library(Matrix) #for bdiag
library(dplyr)
library(tidyverse)
library(pracma)#for pinv
library(xtable)
library("showtext")

## Load fonts
font_add_google("Source Sans Pro", "sourcesanspro")
showtext_auto()


# Load functions
path_func <- "...\\functions"
files <- list.files(path_func, full.names = TRUE)
lapply(files, function(f) {
  message("Sourcing: ", f)
  source(f)
})

## Load data
data <- read_xlsx("TestData1.xlsx")
tau <-0 #Common effects
data$weight <- 1/(data$variance + (tau^2)/2)

### Define variables------------------------------------------------------------
# trials 
trials <- unique(data$trial)
trials <- sort(trials)
M <- length(trials)

# treatments
trts <- unique(data$treatment)
trts <- sort(trts) 
N <- length(trts)

# arms (per trial)
n <- num_arms(data, trials, M)

#sort data so that trials and treatments are always in the same order
data <- data[order(data$treatment), ]
data <- data[order(data$trial), ]


##########################
# Bipartite
##########################

### Bipartite graph-------------------------------------------------------------
# Weighted biadjacency matrix 
B <- bi_adj_matrix(data, M, N, trials, trts)

# Define bipartite graph from adjacency matrix
g.bi <- graph_from_biadjacency_matrix(B, directed=FALSE, weighted=TRUE)

# Plot the graph with layout_as_bipartite
bg <- plot_bi_gg(g.bi, layout = "bipartite", vsize=12, vtxt=20, etxt = 12, endcap=0, startcap=0, emin=0.2, emax=3,
           dodge = 3, push = 12, , border_margin = unit(c(1, 1, 1, 1), "cm"), col.edge = "#707070")

### Bipartite Hat matrix--------------------------------------------------------

## Intermediate matrices
# Within trial covariance matrix
Sig <- get_Sig(data, trials, M, n)
# Between trial covariance matrix
Omega <- get_Omega(data, trials, M, n, tau)
# Weight matrix
W <- solve(Sig+Omega)
# Design matrix 
X <- design_matrix(data, trials, trts, M, N, n)
# Contrast projection matrix
C <- contrast_matrix(data, trials, M, n)

## Hat matrix
# Trial level
H <- solve(t(X)%*%W%*%X)%*%t(X)%*%W
# Arm level
H_arm <- H%*%C

## Consistency equations matrix (to obtain all comparisons)
D <- consistency_matrix(N, trts)

### Bipartite flow network------------------------------------------------------

## Assign flow (H_arm) to each edge in the data
H_arm_full <- D%*%H_arm
for(i in 1:nrow(H_arm_full)){
  cname <- paste("flow", rownames(H_arm_full)[i])
  data[[cname]] <- H_arm_full[i,]
}

## A->B
# Define adjacency matrix for the flow network
FAB <- get_biFlow(H_arm_full, source = "A", sink = "B", N, M, trials, trts, order = 2) #treatment then trial

# Create weighted directed bipartite graph
g.biflow.AB <- graph_from_adjacency_matrix(FAB, mode = "directed", weighted = TRUE, diag = FALSE)
# Set vertex types: FALSE for trials, TRUE for treatments
V(g.biflow.AB)$type <- c(rep(TRUE, N), rep(FALSE, M)) #treatment then trial

#Layout
lay<-layout_with_sugiyama(g.biflow.AB)$layout
rownames(lay)<-c(trts, trials)
colnames(lay)<-c("x","y")
lay["A","x"]<-0
lay["A","y"]<-4.5
lay["B","x"]<-3.5
lay["B","y"]<-4.5
lay["C","x"]<-2.5
lay["C","y"]<-3
lay["D","x"]<-1
lay["D","y"]<-1
lay["1","x"]<-1.75
lay["1","y"]<-7
lay["2","x"]<-1.5
lay["2","y"]<-5
lay["3","x"]<-0.4
lay["3","y"]<-3
lay["4","x"]<-3
lay["4","y"]<-1
lay["5","x"]<-1.75
lay["5","y"]<-2

bf <- plot_biflow_gg(g.biflow.AB, "A", "B", layout = lay, vsize=12, vtxt=20, ashift=7, etxt = 12, col.edge = "#707070",
                     dodge = 3, push = -2, emin=0.2, emax=2.5, border_margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))


## Conjecture 1-----------------------------------------------------------------

## Calculate I' (electrical edge currents)
# Resistance
R <- get_R(data)
# Incidence matrix of bipartite graph
B_bi <- get_bi_incidence(data, trials, trts)
# Nodal currents
Jt <- get_Jt(trts, trials)

Rinv <- solve(R)
L <- t(B_bi)%*%Rinv%*%B_bi

## Edge currents
It <- Jt %*% pinv(L) %*% t(B_bi) %*% Rinv

#####
# Unipartite
#####

## Aggregate hat matrix --------------------------------------------------------
## Intermediate matrices
# Adjusted weight matrix
Wadj <- get_Wadj(data, trials, M, n, tau)
# Aggregate weight matrix
Wa <- get_W_agg(Wadj, trts, N)
# Incidence matrix of unipartite graph
Ba <- get_B_agg(Wa, trts, N)

L <- t(Ba)%*%Wa%*%Ba
Lp <- L_pinv(L)
CN <- get_CN(N)
rownames(CN) <- paste(trts[1], trts[2:N])
colnames(CN) <- trts

## Aggregate Hat matrix
Ha <- CN%*%Lp%*%t(Ba)%*%Wa


### Unipartite graph------------------------------------------------------------

g.uni <- bipartite_projection(g.bi, multiplicity = FALSE, which = "true")
g.uni <- assign_agg_weights(g.uni, Wa)

## Layout (unipartite)
lay.uni<-layout_in_circle(g.uni)
rownames(lay.uni)<-trts
colnames(lay.uni)<-c("x","y")
lay.uni["A","x"]<--1
lay.uni["A","y"]<-1
lay.uni["B","x"]<-1
lay.uni["B","y"]<-1
lay.uni["C","x"]<-1
lay.uni["C","y"]<--1
lay.uni["D","x"]<--1
lay.uni["D","y"]<--1

ug <- plot_uni_gg(g.uni, layout = lay.uni, vsize=12, vtxt=20, etxt = 12, endcap=0, startcap=0, 
            dodge = 3, push = 7, border_margin = unit(c(1, 1, 1, 1), "cm"), col.edge = "#707070",
            emin=0.2, emax=3)

### Unipartite flow network-----------------------------------------------------
Ha_full <- D%*%Ha

FAB.uni <- get_uniFlow(Ha_full, source = "A", sink = "B", N)
g.uniflow.AB <- graph_from_adjacency_matrix(FAB.uni, mode = "directed", weighted = TRUE, diag = FALSE)

uf <- plot_uniflow_gg(g.uniflow.AB, source = "A", sink = "B", layout = lay.uni,
                      vsize=12, vtxt=20, etxt = 12, ashift=7, dodge = 3, push=7, 
                      border_margin = unit(c(1, 1, 1, 1), "cm"), col.edge = "#707070", 
                      emin=0.2, emax=2.5)


########
# Random walks
########

### Random Walks: Conjecture 2--------------------------------------------------

### Transition matrices

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
