library(microbenchmark)
library(Rcpp)

A <- matrix(rnorm(1000000), 100, 10000)
B <- matrix(rnorm(10000), 100, 100)
dim(A); dim(B)

sourceCpp("test.cpp")
microbenchmark(A%*%B,
               armaMatMult(A, B), 
               eigenMatMult(A, B), 
               eigenMapMatMult(A, B), times = 10)
