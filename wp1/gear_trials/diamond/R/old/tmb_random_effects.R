##

library(TMB)
library(ggplot2)
theme_set(theme_bw())

sigma1 <- 2
sigma0 <- 0.1
n <- 20
m <- 10

mu <- rnorm(n, mean = 0, sd = sigma1)

y <- sapply(mu, FUN = function(x){rnorm(m, mean = x, sd = sigma0)})

dat <- data.frame(y = c(y), gp = rep(1:n, each = m))

ggplot(dat, aes(x = gp, y = y)) + geom_point()

## TMB

tmb_model <-"
  #include <TMB.hpp>
  template<class Type>
  Type objective_function<Type>::operator() () {
  // DATA
  DATA_VECTOR(y);
  DATA_VECTOR(gp);
  DATA_SCALAR(m);
  // PARAMETERS:
  PARAMETER(mu);
  PARAMETER_VECTOR(u)
  PARAMETER(log_sigma_obs);
  PARAMETER(log_sigma_proc);
  // PROCEDURE
  Type sigma_proc= exp(log_sigma_proc);
  Type sigma_obs= exp(log_sigma_obs);
  int n = y.size(); 
  Type nll = 0.0;
  // process model:
  for(int i = 0; i < m; i++){
    nll -= dnorm(u[i], Type(0), sigma_proc, true);
  }
  //
  for(int i = 0; i < n; i++){
    Type pred_mean = mu + u(gp(i));
    nll -= dnorm(y[i], pred_mean, sigma_obs, true);
  }
  return nll;
  }
"

write(tmb_model, file = "random_effects.cpp")

compile("random_effects.cpp")

dyn.unload(dynlib("random_effects"))
dyn.load(dynlib("random_effects"))


## SANDBOX


norm.nll <- function(theta, data){
  nll <- - sum(dnorm(data, mean = theta[1], sd = exp(theta[2]), log = TRUE))
}

nsims <- 1e3
sim.df <- data.frame(mean = rep(NA, nsims),
                     tmb.est = NA,
                     r.est = NA)

newtonOption("tol" = 1e-16)

for(i in 1:nsims){
  print(i)
  y <- rnorm(100)
  sim.df$mean[i] <- mean(y)
  ## tmb
  obj <- MakeADFun(
           data = list(y = y),
           parameters = list(mu = 0, log_sigma = 0),
           method = "BFGS",
           DLL = "random_effects",
           inner.control = list(tol10 = 1e-8))
  opt <- do.call("optim", obj)
  sim.df$tmb.est[i] <- opt$par[1]
  rm(obj)
  ## R
  fit <- optim(par = c(0,0), fn = norm.nll, method = "BFGS", data = y)
  sim.df$r.est[i] <- fit$par[1]
}

sqrt(mean((sim.df$tmb.est - sim.df$mean)^2))
sqrt(mean((sim.df$r.est - sim.df$mean)^2))


plot(sim.df[, 1] / sim.df[, 2])

fit <- lm(y ~ 1, data = dat)
coef(fit); opt$par[1]

summary(fit)$sigma; exp(opt$par[2])

