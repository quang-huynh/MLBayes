## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE-------------------------------------------------------
#  browseVignettes("MLZBayes")

## ---- echo = FALSE, fig.cap = 'Figure 1. Traceplot of model parameters from MCMC iterations.', fig.height = 5, fig.width = 7----
library(MLZ)
library(MLZBayes)
library(ggplot2)
library(rstan)
library(ggplot2)
library(loo)
data(Goosefish)
default_priors <- new("MLZ_prior", ncp = 2)

default_priors_1cp <- new("MLZ_prior", ncp = 1)
result.stan <- ML_stan(Goosefish, default_priors)
result.stan.1cp <- ML_stan(Goosefish, default_priors_1cp)
rstan::traceplot(result.stan, pars = c('Z', 'p', 'sigma'))

## ---- echo = FALSE, fig.cap = 'Figure 2. Marginal posterior distributions of model parameters.', fig.height = 5, fig.width = 7----
rstan::stan_hist(result.stan, pars = c('Z', 'p', 'sigma'), bins = 30)

## ---- echo = FALSE, fig.cap = "Figure 3. Joint posterior distribution of the change points in the application to goosefish in the 2-change point model. Lighter blue color indicates a higher probability density.", fig.height = 5, fig.width = 7----
z <- loo::extract_log_lik(result.stan, parameter_name = 'D') + 1963
ggplot(data.frame(z), aes(X1, X2)) + stat_bin2d(bins = 250) + theme_bw() + 
  labs(x = "First change point", y = 'Second change point')

## ---- echo = FALSE, fig.cap = "Figure 4. Observed (black) and posterior predicted (red) mean lengths from the two-change point model (left) and the one-change point model (right).", fig.height = 3, fig.width = 7----
goosefish <- read.csv('more_goosefish.csv')[1:40, ]
par(mfrow = c(1,2), mar = c(5, 4, 1, 1))

res <- rstan::summary(result.stan)$summary

Lpred.ind <- grep('Lpred', rownames(res))
#Lpred.median <- res[Lpred.ind, 6]
#Lpred.2.5 <- res[Lpred.ind, 4]
#Lpred.97.5 <- res[Lpred.ind, 8]
Lpred.mean <- res[Lpred.ind, 1]

plot(goosefish$year, goosefish$mlen, typ = 'o', pch = 16, 
     ylab = "Mean length (cm)", xlab = "Year")
lines(goosefish$year, Lpred.mean, lwd = 3, col = 'red')
#lines(goosefish$year, Lpred.median, lwd = 3, col = 'red')
#lines(goosefish$year, Lpred.2.5, lwd = 3, lty = 2, col = 'red')
#lines(goosefish$year, Lpred.97.5, lwd = 3, lty = 2, col = 'red')

res <- rstan::summary(result.stan.1cp)$summary

Lpred.ind <- grep('Lpred', rownames(res))
#Lpred.median <- res[Lpred.ind, 6]
#Lpred.2.5 <- res[Lpred.ind, 4]
#Lpred.97.5 <- res[Lpred.ind, 8]
Lpred.mean <- res[Lpred.ind, 1]

plot(goosefish$year, goosefish$mlen, typ = 'o', pch = 16, 
     ylab = "Mean length (cm)", xlab = "Year")
lines(goosefish$year, Lpred.mean, lwd = 3, col = 'red')
#lines(goosefish$year, Lpred.median, lwd = 3, col = 'red')
#lines(goosefish$year, Lpred.2.5, lwd = 3, lty = 2, col = 'red')
#lines(goosefish$year, Lpred.97.5, lwd = 3, lty = 2, col = 'red')

