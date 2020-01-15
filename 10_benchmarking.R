library(microbenchmark)
library(Rcpp)

A <- matrix(rnorm(1000000), 100, 10000)
B <- matrix(rnorm(10000), 100, 100)
dim(A); dim(B)

sourceCpp("test.cpp")
TA <- t(A)
TA%*%B


microbenchmark(TA%*%B,
               armaMatMult(TA, B), 
               eigenMatMult(TA, B), 
               eigenMapMatMult(TA, B), times = 10)
