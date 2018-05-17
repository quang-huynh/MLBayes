
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
  
  real<lower=0> muZ;
  real<lower=0> sdZ;
  real<lower=0> musigma;
  real<lower=0> sdsigma; 
}

parameters {
  real<lower=0> Z;
  real<lower=0> sigma;
}


model {
  real Lpred;
  
  // Priors for Z, yearZ, sigma
  Z ~ lognormal(lognormal_mu(muZ, square(sdZ)), lognormal_sd(muZ, square(sdZ))) T[0, 3];
  sigma ~ lognormal(lognormal_mu(musigma, square(sdsigma)), lognormal_sd(musigma, square(sdsigma)));
  
  Lpred = Linf * (1 - (Z/(Z+K)) * (1 - Lc/Linf));
  for(m in 1:count) {
    if(ss[m]>0) target += normal_lpdf(Lobs[m] | Lpred, sigma/sqrt(ss[m]));
  }
}

generated quantities {
  real Lpred;
  vector[count] log_lik;
  Lpred = Linf * (1 - (Z/(Z+K)) * (1 - Lc/Linf));
  for (m in 1:count) log_lik[m] = ss[m] > 0 ? normal_lpdf(Lobs[m] | Lpred, sigma/sqrt(ss[m])) : 0;
}
