library(rstan)

growth.stan <- "
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
"

sink("growth.stan")
cat(growth.stan)
cat("\n")
sink()


# Set up the data 
hake_bayes = hake_df[complete.cases(hake_df[ , c('age', 'length')]), ]
hake_bayes$year_index <- match(hake_bayes$catch_year, unique(hake_bayes$catch_year)) # create a county index with values from 1 to the number of counties in the data set
head(hake_bayes)



hake_bayes_data <- list(y = hake_bayes$length, x = hake_bayes$age, group_index = hake_bayes$year_index, N = nrow(hake_bayes), N_groups = max(hake_bayes$year_index))
subset = sample(1:46596, 300, replace=FALSE)
hake_bayes_data_subset <- list(y = hake_bayes$length[c(subset)], x = hake_bayes$age[c(subset)], group_index = hake_bayes$year_index[c(subset)], N = 300, N_groups = max(hake_bayes$year_index[c(subset)]))


fit <- stan(file = "growth.stan", data = hake_bayes_data, warmup = 1000, iter = 2000)
traceplot(fit2, pars = c("mu_rho", "L1", "sigma_y", "L2", "sigma_rho"), inc_warmup=FALSE, nrow=4)


