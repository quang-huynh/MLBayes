
functions {
#include /functions/misc.stan
}

data {
  int<lower=0> count;
  vector<lower=0>[count] Lobs;
  vector<lower=0>[count] ss;
  real<lower=0> Lc;  
  real<lower=0> Linf;
  real<lower=0> K;
  
  int<lower=0, upper=1> Z_dist;
  vector<lower=0>[2] Z_par;
  int<lower=0, upper=1> sigma_dist;
  vector<lower=0>[2] sigma_par;
}

parameters {
  real<lower=0> Z;
  real<lower=0> sigma;
}

model {
  real Lpred;
  
  // Priors
  if (Z_dist == 0) Z ~ uniform(Z_par[1], Z_par[2]);
  if (Z_dist == 1) Z ~ lognormal(lognormal_mu(Z_par[1], Z_par[2]), lognormal_sd(Z_par[1], Z_par[2])) T[0, 3];
  
  if (sigma_dist == 0) sigma ~ uniform(sigma_par[1], sigma_par[2]);
  if (sigma_dist == 1) sigma ~ lognormal(lognormal_mu(sigma_par[1], sigma_par[2]), lognormal_sd(sigma_par[1], sigma_par[2]));
  
  // Generate predicted mean length 
  Lpred = Linf * (1 - (Z/(Z+K)) * (1 - Lc/Linf));
  
  // Likelihood
  for (m in 1:count) {
    if (ss[m]>0) target += normal_lpdf(Lobs[m] | Lpred, sigma/sqrt(ss[m]));
  }
}

generated quantities {
  real Lpred;
  vector[count] log_lik;
  Lpred = Linf * (1 - (Z/(Z+K)) * (1 - Lc/Linf));
  for (m in 1:count) log_lik[m] = ss[m] > 0 ? normal_lpdf(Lobs[m] | Lpred, sigma/sqrt(ss[m])) : 0;
}
