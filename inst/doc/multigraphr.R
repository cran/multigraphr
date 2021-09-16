## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(out.width = "100%",
  cache = FALSE
)

## ----eval=TRUE----------------------------------------------------------------
library('multigraphr')

## ----adj, include=TRUE, results='markup', message=FALSE-----------------------
A <-  matrix(c(1, 1, 0, 
               1, 2, 2, 
               0, 2, 0), 
             nrow = 3, ncol = 3)
A

## ----degseq, include=TRUE, results='markup', message=FALSE--------------------
D <- get_degree_seq(adj = A, type = 'graph')
D

## ----rsm_ex1, include=TRUE, results='markup', message=FALSE-------------------
rsm_1 <- rsm_model(deg.seq = D)
rsm_1$m.seq

## ----rsm_ex2, include=TRUE, results='markup', message=FALSE-------------------
rsm_1$prob.dists

## ----ieas_ex1, include=TRUE, results='markup', message=FALSE------------------
ieas_1 <-   iea_model(adj = A , type = 'graph',  model = 'IEAS', K = 0, apx = TRUE)
isa_1 <-   iea_model(adj = A , type = 'graph',  model = 'ISA', K = 0, apx = TRUE)
isa_1$nr.multigraphs
ieas_1$nr.multigraphs

## ----ieas_ex3, include=TRUE, results='markup', message=FALSE------------------
ieas_2 <-   iea_model(adj = A , type = 'graph', model = 'IEAS', 
                      K = 0, apx = FALSE)

## ----isa_ex2, include=TRUE, results='markup', message=FALSE-------------------
isa_2 <-   iea_model(adj = A , type = 'graph', model = 'ISA', 
                     K = 0, apx = FALSE, p.seq = c(1/3, 1/3, 1/3))

## ----rsm_ex3, include=TRUE, results='markup', message=FALSE-------------------
rsm_1$M

## ----ieas_ex2, include=TRUE, results='markup', message=FALSE------------------
ieas_1$M
ieas_1$R

## ----isa_ex1, include=TRUE, results='markup', message=FALSE-------------------
isa_1$M
isa_1$R

## ----gof1, include=TRUE, results='markup', message=FALSE, eval=FALSE----------
#  gof1 <- gof_sim(m = 10, model = 'IEAS', deg.mod = c(8,8,2,2),
#                  hyp = 'IEAS', deg.hyp = c(6,6,6,2))

## ----gof2, include=TRUE, results='markup', message=FALSE, eval=FALSE----------
#  gof2 <- gof_sim(m = 10, model = 'IEAS', deg.mod = c(14,2,2,2),
#                  hyp = 'IEAS', deg.hyp = c(14,2,2,2))

## ----gof3, include=TRUE, results='markup', message=FALSE, eval=FALSE----------
#  gof3 <- gof_sim(m = 10, model = 'RSM', deg.mod = c(14,2,2,2),
#                  hyp = 'IEAS', deg.hyp = 0)

## ----gof4, include=TRUE, results='markup', message=FALSE, eval=FALSE----------
#  gof4 <- gof_sim(m = 10, model = 'ISA', deg.mod = c(14,2,2,2),
#                  hyp = 'ISA', deg.hyp = 0)

## ----flor1, include=TRUE, results='markup', message=FALSE---------------------
flor_m <- t(matrix(c (0, 0, 1, 0, 0, 0,	0, 0,
                      0, 0, 0, 0, 0, 0,	0, 0,
                      0, 0,	0, 2, 0, 0,	1, 5,
                      0, 0,	0, 0, 0, 0,	1, 1,
                      0, 0,	0, 0, 0, 0,	1, 2,
                      0, 0,	0, 0, 0, 0,	2, 1,
                      0, 0,	0, 0, 0, 0,	0, 2,
                      0, 0,	0, 0, 0, 0,	0, 1), nrow= 8, ncol=8))

## ----flor2, include=TRUE, results='markup', message=FALSE---------------------
flor_adj <- flor_m+t(flor_m)
flor_adj 

## ----flor3, include=TRUE, results='markup', message=FALSE---------------------
flor_d <- get_degree_seq(adj = flor_adj, type = 'multigraph')
flor_d

## ----flor4, include=TRUE, results='markup', message=FALSE---------------------
flor_ieas_test <- gof_test(flor_adj, 'multigraph', 'IEAS', flor_d, 35)
flor_ieas_test

## ----flor5, include=TRUE, results='markup', message=FALSE---------------------
flor_isa_test <- gof_test(flor_adj, 'multigraph', 'ISA', flor_d, 35)
flor_isa_test 

## ----func1, include=TRUE, results='markup', message=FALSE---------------------
r <- (2*3)/2 # vertex pair sites (or length of edge multiplicity sequences)
mg <- nsumk(r,4) # number of rows give number of possible multigraphs
mg

## ----func2, include=TRUE, results='markup', message=FALSE---------------------
Q <- get_edge_assignment_probs(m = 8, deg.seq = c(4,4,4,4), model = 'IEAS')
Q

## ----func3, include=TRUE, results='markup', message=FALSE---------------------
mg <- get_edge_multip_seq(deg.seq = c(4,2,2))
mg

