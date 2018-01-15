library(rstan)

## create data set for 2cmt example & visualize results

options(rcpp.cache.dir=".")

expose_stan_functions("oral_2cmt_mpi4.stan")

dose <- 5
lstate0 <- c(log(dose), -35, -35)

Dt1 <- round(exp(seq(log(0.5), log(5), length=6)), 2)
Dt2 <- round(exp(seq(log(7), log(24*7), length=11)), 2)
Dt <- c(Dt1, Dt2)


## k = CL/V
## MRT = 1/k = V/CL
## V = MRT * CL

## sqrt(V1*V2) = V
lV <- log(8)
## sqrt(CL1*CL2) = CL
lCL <- lV + log(log(2)/10)

## sqrt(V1/V2)
lVr <- 0.5 * log(0.2)

## sqrt(CL/Q)
lCLr <- 0.5 * log(0.3)

lV1 <- lV + lVr
lV2 <- lV - lVr

lCL1 <- lCL + lCLr
lCL2  <- lCL - lCLr

lka <- log(log(2)/2)

lfrel <- 0

## note: the Stan program fixes the relative bio-availability to 1 for
## the population

##true_theta <- c(lfrel, lka, lCL, lCLr, lV, lVr)
true_theta <- c(lfrel, lka, lCL, lCLr, lV, lVr)

true_omega <- c(0, 0.1, 0.3, 0, 0.2, 0)

pk_sys <- function(lstate0, theta) {
    function(x) {
        twoCmtOral_analytic(lstate0[1] + theta[1], x, theta[-1])
    }
}

sys <- pk_sys(lstate0, true_theta)

curve(exp(sys(x)), 0.01, max(Dt))
points(Dt, exp(sys(Dt)))

curve(sys(x), 0.1, max(Dt))
points(Dt, sys(Dt))

## creates a stan data set
make_ds <- function() {
    stan_data_base <- list(
                      T=length(Dt),
                      time=Dt,
                      true_theta=true_theta,
                      true_omega=true_omega,
                      use_ode=0,
                      known_omega=1,
                      known_sigma_y=0,
                      dose=dose
                      )

    stan_rdump(names(stan_data_base),
               paste0("oral4_stan-base.R"),
               envir=list2env(stan_data_base))
}

make_ds()

## visualize results
library(bayesplot)

J <- 50

runs <- list(mpi=read_stan_csv(paste0("samples5-1-rect-", J, "-0.csv")),
             serial=read_stan_csv(paste0("samples5-1-rect-", J, "-1.csv")),
             stan=read_stan_csv(paste0("samples5-1-rect-", J, "-2.csv")))

runs <- list(mpi=read_stan_csv(paste0("samples5-1-rect-", J, "-0.csv")))

total <- colSums(sapply(runs, get_elapsed_time))

## runtime in minutes
round(total/60, 2)
## speedup vs vanilla stan version
round(total["stan"]/total, 2)
## speedup vs serial version
round(total["serial"]/total, 2)

bias <- lapply(runs, extract, pars=c("bias_theta", "bias_sigma_y"), inc_warmup=FALSE, permuted=FALSE)

flat <- list()
for(i in names(runs)) {
    r <- as.matrix(bias[[i]][,1,])
    colnames(r) <- paste(i, colnames(r), sep="/")
    flat <- c(list(r), flat)
}

flat <- do.call(cbind, flat)

batches <- gsub("[a-z]+\\/", "", colnames(flat))
version <- gsub("\\/.+", "", colnames(flat))

version
colnames(flat) <- version

mcmc_recover_intervals(flat, rep(0, ncol(flat)), batch=batches) + facet_text() + hline_0()

print(runs[[1]], pars=c("theta_v", "sigma_y"))

sr <- runs[[1]]

traceplot(sr, pars="theta_raw", inc_warmup=TRUE)

traceplot(sr, pars="theta_raw", inc_warmup=FALSE)

pairs(sr, pars="theta_raw")

sp <- get_sampler_params(sr)[[1]]

Nwarmup <- nrow(sp)/2
dim(sp)

spi <- sp[- (1:Nwarmup),]

head(sp)
colSums(spi)

colMeans(spi)

spi[spi[,"divergent__"]==1,]
