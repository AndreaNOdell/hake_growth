
data {
  int<lower=0> N;   // number of data items
  int<lower=0> N_groups; // number of groups
  vector[N] x;   // predictor vector (age)
  vector[N] y;      // outcome vector (length)
  int<lower=1,upper=18> group_index[N]; // group index
}
parameters {
  real rho[N_groups];	// random K
  real mu_rho;        // overall intercept
  real L1;            // size at age 2
  real L2;            // size at age 15
  real<lower=0> sigma_y;  // error scale -- individual level
  real<lower=0> sigma_rho;  // error scale -- group level
}
model {
  mu_rho ~ lognormal(0, 4); // priors
  L1 ~ normal(32,5);
  L2 ~ normal(51,5);
  sigma_y ~ gamma(2, 0.5);
  sigma_rho ~ gamma(2, 0.5);
  for (n in 1:N_groups)
	  rho[n] ~ normal(mu_rho, sigma_rho); // likelihood, group level
  for (n in 1:N)
  	y[n] ~ normal(L1 + (L2-L1)*((1-(rho[group_index[n]]^(x[n]-2)))/(1-(rho[group_index[n]]^(13)))), sigma_y);  // likelihood, individual level
}

