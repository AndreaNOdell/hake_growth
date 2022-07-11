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


growth.vb.stan <- "
data {
  int<lower=0> N;   // number of data items
  int<lower=0> N_groups; // number of groups
  vector[N] x;   // predictor vector (age)
  vector[N] y;      // outcome vector (length)
  int<lower=1,upper=N_groups> group_index[N]; // group index
}
parameters {
  vector[N_groups] K;	// random K
  real mu_K;        // overall intercept
  real Linf;            // size at max age
  real L0;            // size at age 0
  real<lower=0> sigma_y;  // error scale -- individual level
  real<lower=0> sigma_K;  // error scale -- group level
}
model {
  mu_K ~ uniform(0, 4); // priors
  Linf ~ normal(50,5);
  L0 ~ normal(1,0.001);
  sigma_y ~ gamma(2, 0.5);
  sigma_K ~ gamma(2, 0.5);
  for (n in 1:N_groups)
	  K[n] ~ normal(mu_K, sigma_K); // likelihood, group level
  for (i in 1:N)
  	y[i] ~ normal(Linf- (Linf - L0) * exp(-K[group_index[i]] * x[i]) , sigma_y);  // likelihood, individual level
}
"

sink("growth.vb.stan")
cat(growth.vb.stan)
cat("\n")
sink()



# Set up the data 
hake_bayes = hake_df[complete.cases(hake_df[ , c('age', 'length')]), ]
hake_bayes$year_index <- match(hake_bayes$catch_year, unique(hake_bayes$catch_year)) # create a county index with values from 1 to the number of counties in the data set
head(hake_bayes)



hake_bayes_data <- list(y = hake_bayes$length, x = hake_bayes$age, group_index = hake_bayes$year_index, N = nrow(hake_bayes), N_groups = max(hake_bayes$year_index))
subset = sample(1:46596, 300, replace=FALSE)
hake_bayes_data_subset <- list(y = hake_bayes$length[c(subset)], x = hake_bayes$age[c(subset)], group_index = hake_bayes$year_index[c(subset)], N = 300, N_groups = max(hake_bayes$year_index[c(subset)]))


stanfit <- stan(file = "growth.vb.stan", data = hake_bayes_data, warmup = 1000, iter = 5000)
#traceplot(fit2, pars = c("mu_rho", "L1", "sigma_y", "L2", "sigma_rho"), inc_warmup=FALSE, nrow=4)




## That didn't work well, but I found a package called BayesGrowth that will make fitting
# Bayesian growth models simpler... so let's give it a try. The caveat is that it doesn't allow
# hierarchical/mixed effects

library(BayesGrowth)

hake_bayes = hake_df[complete.cases(hake_df[ , c('age', 'length')]), ] %>% 
  select(age, length) %>% 
  rename(Age = age, Length = length)

## Biological info - lengths in mm
max_size <- 60
max_size_se <- 5
birth_size <- 1
birth_size_se <- 0.001 # an se cannot be zero

# Use the function to estimate the rstan model
fit <- Estimate_MCMC_Growth(data = hake_bayes, 
                            Model = "VB" ,
                            iter = 5000,
                            Linf = max_size,
                            Linf.se = max_size_se,
                            L0 = birth_size,
                            sigma.max = 100,
                            L0.se = birth_size_se,
                            k.max = 1)

list_of_draws <- extract(fit,c("Linf", "k","L0", "sigma")) %>% 
  as.data.frame() %>% 
  gather(Parameter, Value) %>% 
  filter(Parameter %in% c("Linf", "k","L0", "sigma"))


ggplot(list_of_draws, aes(Value))+
  geom_density(fill = "royalblue")+
  facet_wrap(~Parameter, scales = "free", ncol = 2)+
  theme_bw()








