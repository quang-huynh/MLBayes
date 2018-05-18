
  real lognormal_mu(real mean, real sd) {
    return log(mean/(sqrt(1 + sd*sd/(mean * mean))));
  }
  
  real lognormal_sd(real mean, real sd) {
    return sqrt(log(1 + sd*sd/(mean * mean)));
  }
  