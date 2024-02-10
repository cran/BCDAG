## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
oldpar <- par(no.readonly = TRUE)
oldoptions <- options()

## ----setup--------------------------------------------------------------------
library(BCDAG)

## -----------------------------------------------------------------------------
## Generate data
set.seed(1)
q <- 8
w <- 0.2
DAG <- rDAG(q,w)
a <- q
U <- diag(1,q)
outDL <- rDAGWishart(n=1, DAG, a, U)
L <- outDL$L; D <- outDL$D
Omega <- L %*% solve(D) %*% t(L)
Sigma <- solve(Omega)
n <- 1000
X <- mvtnorm::rmvnorm(n = n, sigma = Sigma)

## ----echo = FALSE, include=FALSE----------------------------------------------
out <- learn_DAG(S = 5000, burn = 1000, data = X,
                 a, U, w, 
                 fast = FALSE, save.memory = FALSE, collapse = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  ## Run MCMC
#  out <- learn_DAG(S = 5000, burn = 1000, data = X,
#                   a, U, w,
#                   fast = FALSE, save.memory = FALSE, collapse = FALSE)

## ----fig.width = 7, fig.height= 6---------------------------------------------
get_diagnostics(out, ask = FALSE)

## -----------------------------------------------------------------------------
print(out)
summary(out)
plot(out)

## -----------------------------------------------------------------------------
get_edgeprobs(out)

## -----------------------------------------------------------------------------
MPMdag <- get_MPMdag(out)
MPMdag

## -----------------------------------------------------------------------------
MAPdag <- get_MAPdag(out)
MAPdag

## ----fig.width = 7------------------------------------------------------------
par(mfrow = c(1,3))
Rgraphviz::plot(as_graphNEL(DAG), main = "True DAG")
Rgraphviz::plot(as_graphNEL(MPMdag), main = "MPM DAG")
Rgraphviz::plot(as_graphNEL(MAPdag), main = "MAP DAG")

## -----------------------------------------------------------------------------
round(L, 3)

## -----------------------------------------------------------------------------
Rgraphviz::plot(as_graphNEL(DAG), main = "True DAG")

## -----------------------------------------------------------------------------
causaleffect(targets = c(4,5), response = 1, L = L, D = D)
causaleffect(targets = 4, response = 1, L = L, D = D)
causaleffect(targets = 5, response = 1, L = L, D = D)

## -----------------------------------------------------------------------------
DAG2 <- DAG
DAG2[4,5] <- 1
par(mfrow = c(1,2))
Rgraphviz::plot(as_graphNEL(DAG), main = "True DAG")
Rgraphviz::plot(as_graphNEL(DAG2), main = "Modified DAG")

## ----include = FALSE----------------------------------------------------------
par(mfrow = c(1,1))

## -----------------------------------------------------------------------------
L2 <- L
L2[4,5] <- runif(1)
L2[4,5]

## -----------------------------------------------------------------------------
causaleffect(targets = c(4,5), response = 1, L = L2, D = D)
causaleffect(targets = 4, response = 1, L = L2, D = D)
causaleffect(targets = 5, response = 1, L = L2, D = D)

## -----------------------------------------------------------------------------
effects_out <- get_causaleffect(out, targets = c(4,5), response = 1)
head(effects_out$causaleffects)

## -----------------------------------------------------------------------------
print(effects_out)
summary(effects_out)
plot(effects_out)

## ----echo = FALSE, include=FALSE----------------------------------------------
coll_out <- learn_DAG(S = 5000, burn = 1000, data = X,
                      a, U, w, 
                 fast = FALSE, save.memory = FALSE, collapse = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  coll_out <- learn_DAG(S = 5000, burn = 1000, data = X,
#                        a, U, w,
#                   fast = FALSE, save.memory = FALSE, collapse = TRUE)

## -----------------------------------------------------------------------------
names(coll_out)

## ----include = FALSE----------------------------------------------------------
effects_collout <- get_causaleffect(coll_out, targets = c(4,5), response = 1)

## ----eval = FALSE-------------------------------------------------------------
#  effects_collout <- get_causaleffect(coll_out, targets = c(4,5), response = 1)

## -----------------------------------------------------------------------------
round(effects_collout$post_mean, 3)

## ----include = FALSE----------------------------------------------------------
par(oldpar)
options(oldoptions)

