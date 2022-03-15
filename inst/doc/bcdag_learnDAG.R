## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
oldpar <- par(no.readonly = TRUE)
oldoptions <- options()

## ----setup--------------------------------------------------------------------
library(BCDAG)

## -----------------------------------------------------------------------------
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
#  out <- learn_DAG(S = 5000, burn = 1000, data = X,
#                   a, U, w,
#                   fast = FALSE, save.memory = FALSE, collapse = FALSE)

## -----------------------------------------------------------------------------
class(out)

## -----------------------------------------------------------------------------
str(out)

## -----------------------------------------------------------------------------
out$Graphs[,,1]
round(out$L[,,1],2)
round(out$D[,,1],2)

## ----echo = FALSE, include=FALSE----------------------------------------------
collapsed_out <- learn_DAG(S = 5000, burn = 1000, data = X,
                 a, U, w, 
                 fast = FALSE, save.memory = FALSE, collapse = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  collapsed_out <- learn_DAG(S = 5000, burn = 1000, data = X,
#                   a, U, w,
#                   fast = FALSE, save.memory = FALSE, collapse = TRUE)

## -----------------------------------------------------------------------------
names(collapsed_out)
class(collapsed_out)
attributes(collapsed_out)$type
collapsed_out$Graphs[,,1]

## ----echo = FALSE, include=FALSE----------------------------------------------
compressed_out <- learn_DAG(S = 5000, burn = 1000, data = X,
                 a, U, w, 
                 fast = FALSE, save.memory = TRUE, collapse = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  compressed_out <- learn_DAG(S = 5000, burn = 1000, data = X,
#                   a, U, w,
#                   fast = FALSE, save.memory = TRUE, collapse = FALSE)

## -----------------------------------------------------------------------------
names(compressed_out)
class(compressed_out)
attributes(compressed_out)$type

## -----------------------------------------------------------------------------
compressed_out$Graphs[1]
compressed_out$L[1]
compressed_out$D[1]

## -----------------------------------------------------------------------------
bd_decode(compressed_out$Graphs[1])
round(bd_decode(compressed_out$L[1]),2)
round(bd_decode(compressed_out$D[1]),2)

## ----echo = FALSE, include=FALSE----------------------------------------------
comprcoll_out <- learn_DAG(S = 5000, burn = 1000, data = X,
                 a, U, w, 
                 fast = FALSE, save.memory = TRUE, collapse = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  comprcoll_out <- learn_DAG(S = 5000, burn = 1000, data = X,
#                   a, U, w,
#                   fast = FALSE, save.memory = TRUE, collapse = TRUE)

## -----------------------------------------------------------------------------
names(comprcoll_out)
class(comprcoll_out)
attributes(comprcoll_out)$type
bd_decode(comprcoll_out$Graphs[1])

## ----results='hide'-----------------------------------------------------------
# No approximation
time_nofast <- system.time(out_nofast <- learn_DAG(S = 5000, burn = 1000, data = X, 
                      a, U, w, 
                      fast = FALSE, save.memory = FALSE, collapse = FALSE))
# Approximation
time_fast <- system.time(out_fast <- learn_DAG(S = 5000, burn = 1000, data = X, 
                      a, U, w, 
                      fast = TRUE, save.memory = FALSE, collapse = FALSE))

## -----------------------------------------------------------------------------
time_nofast
time_fast

## -----------------------------------------------------------------------------
round(get_edgeprobs(out_nofast), 2)
round(get_edgeprobs(out_fast), 2)

## ---- include = FALSE---------------------------------------------------------
par(oldpar)
options(oldoptions)

