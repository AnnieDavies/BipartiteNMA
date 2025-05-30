## Get results for fictional example using netmeta functions
## Gerta Ruecker 2025

library(netmeta)
library(Matrix)

# Read in data
data <- read_xlsx("TestData1.xlsx")

p1 <- pairwise(treat = data$treatment, TE = data$mean, seTE = sqrt(data$variance), studlab = data$trial)
net1 <- netmeta(p1)
netgraph(net1, multiarm = TRUE, cex = 2)
summary(net1)
forest(net1, ref = "A")

# Aggregate hat matrix (eq 18)
Hagg <- round(hatmatrix(net1, method = "Davies")$common[1:3,],3)

# Trial level hat matrix
# Take the Krahn matrix (which must be permuted):
H <- hatmatrix(net1, method = "Krahn", type = "studies")$common[1:3, c(1,4,5,2,6,7,3)]

# Define C:
C2 <- rbind(c(-1, 1))
C3 <- rbind(c(-1, 1, 0),
            c(-1, 0, 1))
C <- bdiag(C2, C3, C2, C3, C2)

# Arm-level hat matrix, H⊥⊤:
Harm <- round(H %*% C,3)

## Check NMA estimates
H %*% C %*% data$mean                      # compare to forest plot:
forest(net1, ref = "A", digits = 7) # or to
net1$TE.nma.common
