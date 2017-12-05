library(rstan)

set.seed(1)
run0 <- read_stan_csv("samples-1-rect-0.csv")
set.seed(1)
run1 <- read_stan_csv("samples-1-rect-1.csv")
set.seed(1)
run2 <- read_stan_csv("samples-1-rect-2.csv")


print(run0, "lp__")
print(run1, "lp__")
print(run2, "lp__")

lp0 <- extract(run0, "lp__", inc_warmup=TRUE, permuted=FALSE)
lp1 <- extract(run1, "lp__", inc_warmup=TRUE, permuted=FALSE)
lp2 <- extract(run2, "lp__", inc_warmup=TRUE, permuted=FALSE)

lp0 - lp1 > 1E-14
lp0 - lp2

sum(abs(lp0 - lp1))
sum(abs(lp0 - lp2))
sum(abs(lp1 - lp2))
