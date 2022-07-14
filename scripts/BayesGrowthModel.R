library(rstan)
options(mc.cores = parallel::detectCores() - 1)
rstan_options(auto_write = TRUE)

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
  vector[N] x;   // predictor vector (age)
  vector[N] y;      // outcome vector (length)
}
parameters {
  real K;	          // growth K
  real Linf;     // random size at max age
  real L0;          // size at age 0
  real<lower=0> sigma_y;  // error scale -- individual level
}
model {
  Linf ~ normal(50,1); // priors
  K ~ uniform(0.3, 0.5); 
  L0 ~ normal(7,1);
  sigma_y ~ normal(3.5, 0.5);
  for (i in 1:N)
  	y[i] ~ normal(Linf - (Linf - L0) * exp(-K * x[i]) , sigma_y);  // likelihood, individual level
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

hake_bayes_simple  <- list(y = hake_bayes$length, x = hake_bayes$age, N = nrow(hake_bayes))

stanfit <- stan(file = "growth.vb.stan", data = hake_bayes_simple, warmup = 1000, iter = 5000)
#traceplot(fit2, pars = c("mu_rho", "L1", "sigma_y", "L2", "sigma_rho"), inc_warmup=FALSE, nrow=4)




## That didn't work well, but I found a package called BayesGrowth that will make fitting
# Bayesian growth models simpler... so let's give it a try. The caveat is that it doesn't allow
# hierarchical/mixed effects

library(BayesGrowth)

hake_bayes = hake_df[complete.cases(hake_df[ , c('age', 'length')]), ] 
hake_bayes$year_index <- match(hake_bayes$catch_year, unique(hake_bayes$catch_year))
hake_bayes = hake_bayes %>% 
  select(age, length) %>% 
  rename(Age = age, Length = length)
head(hake_bayes)

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

# Use the function to estimate the rstan model
fit <- Estimate_MCMC_Growth_mixed(data = hake_bayes, 
                            Model = "VB" ,
                            iter = 5000,
                            Linf = max_size,
                            Linf.se = max_size_se,
                            L0 = birth_size,
                            sigma.max = 100,
                            L0.se = birth_size_se,
                            k.max = 1)

list_of_draws <- extract(fit, c("Linf", "k","L0", "sigma")) %>% 
  as.data.frame() %>% 
  gather(Parameter, Value) %>% 
  filter(Parameter %in% c("Linf", "k","L0", "sigma"))


ggplot(list_of_draws, aes(Value))+
  geom_density(fill = "royalblue")+
  facet_wrap(~Parameter, scales = "free", ncol = 2)+
  theme_bw()

#save(fit, file = "results/BayesGrowthFit.RData")

plot(hake_df$age, hake_df$length, type = "p")
lines(seq(0,25,by = 0.1), 49.71 - (49.71 - 1)*exp(-0.51*seq(0,25,by = 0.1)), col = 'red')

# Linf = 49.71
# k = 0.51
# L0 = 1
# sigma = 3.69



growth.smart.stan <- "
data {
  int<lower=1> n; //number of samples
  //data vectors
  vector<lower=0>[n] age; //age data
  vector<lower=0>[n] length; //length data
}

parameters {
  //VBGM parameters
  real<lower=0> L0; //length-at-birth
  real<lower=0> Linf; //asymptotic length
  real<lower=0> K; // growth coefficient
  
  //Likelihood parameters
  real<lower=0> sigma; //RSE
}

model {
  //storage
  vector[n] PredL; //predicted lengths
  
  //VBGM priors
  Linf ~ normal(51,1);
  L0 ~ normal(22,1);
  K ~ uniform(0.37, 0.4);
  sigma ~ lognormal(0, 2);
  
  //VBGM likelihood
  for(i in 1:n) {
    PredL[i] = Linf - (Linf - L0)*exp(-K*age[i]);
    target += normal_lpdf(length[i]|PredL[i], sigma);//likelihood
    }
}
"

sink("growth.smart.stan")
cat(growth.smart.stan)
cat("\n")
sink()


mean.age<-tapply(hake_bayes$Length, round(hake_bayes$Age), mean,na.rm = T)
Lt1<-mean.age[2:length(mean.age)]
Lt<-mean.age[1:length(mean.age)-1]
model<-lm(Lt1 ~ Lt)
k <- suppressWarnings(abs(-log(model$coef[2]))) #in case a nan occurs
k <- ifelse(is.nan(k),0.1,k) # in case a nan occurs
Linf<-abs(model$coef[1]/(1-model$coef[2]))
L0 <-lm(mean.age ~ poly(as.numeric(names(mean.age)), 2, raw = TRUE))$coef[1]

startingparameters = list(list(Linf = Linf), list(L0 = L0), list(k = k), list(sigma = 3))

hake_bayes = hake_df[complete.cases(hake_df[ , c('age', 'length')]), ]
bayes_data = list(length = hake_bayes$length, age = hake_bayes$age, n = nrow(hake_bayes))

stanfit <- stan("growth.smart.stan", data = bayes_data, warmup = 5000, iter = 10000, init = startingparameters)


list_of_draws <- extract(stanfit,c("Linf", "K","L0", "sigma")) %>% 
  as.data.frame() %>% 
  gather(Parameter, Value) %>% 
  filter(Parameter %in% c("Linf", "K","L0", "sigma"))


ggplot(list_of_draws, aes(Value))+
  geom_density(fill = "royalblue")+
  facet_wrap(~Parameter, scales = "free", ncol = 2)+
  theme_bw()




# Let's try brms.... 
library(brms)

hake_brms = hake_df[complete.cases(hake_df[ , c('age', 'length')]), ] %>% 
  select(age, length, catch_year)

fit <- brm(
  bf(length ~ L1 + (L2 - L1)*((1 - (rho)^(age-2))/(1 - (rho)^(13))),
     L1 ~ 1, L2 ~ 1, rho ~ 1 + (1|catch_year),
     nl = TRUE),
  data = hake_brms, family = gaussian(),
  prior = c(
    prior(normal(32, 3), nlpar = "L1"),
    prior(normal(51, 3), nlpar = "L2"),
    prior(lognormal(0, 2), nlpar = "rho")
  ),
  control = list(adapt_delta = 0.9)
)

# saveRDS(fit, file = "bayes_results/brms_fixed_fit.RData")
load(file = "bayes_results/brms_fixed_fit.RData")

traceplot(fit, pars = c("L1", "L2", "rho"))




