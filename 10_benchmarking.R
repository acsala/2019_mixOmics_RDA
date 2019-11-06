

library(microbenchmark)
mb = microbenchmark(transpose(X), t(X), times = 10, unit = "s")
print(mb)

library(Rcpp)

A <- matrix(rnorm(10000), 100, 100)
B <- matrix(rnorm(10000), 100, 100)

library(microbenchmark)
sourceCpp("test.cpp")
microbenchmark(A%*%B, armaMatMult(A, B), eigenMatMult(A, B), eigenMapMatMult(A, B), times = 10)
