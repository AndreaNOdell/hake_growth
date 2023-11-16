library(tidyverse)
library(brms)
library(mgcv)
library(ggpubr)
library(performance)
library(tidybayes)
library(modelr)
library(bayesplot)
library(gratia)

df = read.csv("raw_data/selection.csv") # load in data
#colnames(df)

hake_df = df %>% # make new dataframe with relevant information only
  dplyr::filter(scientific_name == "Merluccius productus") %>% 
  dplyr::select("age", "avg_weight", "distance_fished", #"hb_date", 
         "length", "sex_description", "weight", "catch_modified_date",
         "eq_date", "hb_date", 'hb_latitude', 'hb_longitude', "haul_num", "duration") %>% 
  tidyr::separate(col = eq_date , into = c('date', 'time'), sep = " ", remove = FALSE) %>% 
  tidyr::separate(col = date, into = c("catch_year", "catch_month", 'catch_day'), sep = "-")

# Set missing catch_year values to the year inputed for the other date related columns
hake_df[is.na(hake_df$catch_month),]$catch_year <- 2017

# create plus group
hake_df = hake_df %>% 
  mutate(new_age = age)
hake_df$new_age = as.integer(lapply(hake_df$new_age, function(x) ifelse(x > 14, 15, x)))
hake_df$cohort = as.numeric(hake_df$catch_year) - hake_df$new_age

ggplot(hake_df, aes(x = new_age, y = weight)) +
  geom_point() +
  theme_classic()


# Create a dataframe with age and weight, complete cases only
hake_weight_age_df = hake_df %>% 
  dplyr::select(new_age, weight, catch_year, catch_month, catch_day, cohort, length, distance_fished, sex_description, hb_latitude, hb_longitude) %>% 
  tidyr::drop_na()
hake_weight_age_df$catch_month = as.numeric(hake_weight_age_df$catch_month)
hake_weight_age_df$catch_day = as.numeric(hake_weight_age_df$catch_day)
#save(hake_weight_age_df, file = "results/hake_weight_age_df.RData")


intercept_only_brms_out <- brm(bf(weight ~ 1 + (1 | new_age)),
              data = hake_weight_age_df,                               
              warmup = 1000, 
              iter   = 3000, 
              chains = 4, 
              inits  = "random",
              cores  = 2)
estimated_intercepts = coef(intercept_only_brms_out)$new_age[,1,]
#save(intercept_only_brms_out, file = "results/intercept_only_brms_out.RData")

plot(estimated_intercepts, xlab = "age", ylab = "estimated intercepts")


# gam ----------------------------------------------------
gam_out <- gam(weight ~ s(new_age), data = hake_weight_age_df, method = "REML")
#save(gam_out, file = "results/gam_out.RData")
plot(gam_out, xlab = "age")

# make catch_year a factor
hake_weight_age_df$catch_year = as.factor(hake_weight_age_df$catch_year)
hake_weight_age_df$cohort = as.factor(hake_weight_age_df$cohort)

# group level smoothers with different wiggliness
gamm_out = gam(weight ~ s(new_age, bs="tp") +
      s(new_age, by = catch_year, m = 1, bs="tp") +
      s(catch_year, bs="re"), 
      data = hake_weight_age_df, method="REML")


# group level smoothers with same wiggliness
gamm_GS_out = gam(weight ~ s(new_age, m = 2) +
                 s(new_age, catch_year, bs="fs", m = 2), 
               data = hake_weight_age_df, method="REML")

gamm_year_out = gam(weight ~ s(new_age) + s(catch_year, bs="re"), 
                  data = hake_weight_age_df, method="REML")
gamm_year_cohort_out = gam(weight ~ s(new_age) + s(catch_year, bs="re")
                           + s(cohort, bs="re"), data = hake_weight_age_df, method="REML")
gamm_cohort_out = gam(weight ~ s(new_age) + s(cohort, bs="re"), data = hake_weight_age_df, method="REML")


gamm_year_cohort_spatial_out = gam(weight ~ s(new_age) + s(catch_year, bs="re") 
                                   + te(hb_longitude,hb_latitude,bs="tp")
                                   + s(cohort, bs="re"), data = hake_weight_age_df, method="REML")

gamm_trial_out = gamm(weight ~ s(new_age) + s(catch_year, bs="re")
                           + s(new_age, cohort, bs="fs"), data = hake_weight_age_df, method="REML")

# for cohorts
# group level smoothers with different wiggliness
gamm_cohort_out = gam(weight ~ s(new_age, bs="tp") +
                        s(new_age, by = cohort, m = 1, bs="tp") +
                        s(cohort, bs="re"), 
                      data = hake_weight_age_df, method="REML")


# group level smoothers with same wiggliness
gamm_GS_cohort_out = gam(weight ~ s(new_age, m = 2) +
                    s(new_age, cohort, bs="fs", m = 2), 
                  data = hake_weight_age_df, method="REML")


# group level smoothers with same wiggliness for both year and cohort
gamm_GS_year_cohort_out = gam(weight ~ s(new_age, m = 2) +
                           s(new_age, cohort, bs="fs", m = 2) +
                           s(new_age, catch_year, bs="fs", m = 2) , 
                         data = hake_weight_age_df, method="REML")

# save residuals into dataframe
hake_weight_age_df$gam.resids = residuals(gam_out, type = "deviance")
hake_weight_age_df$gamm.resids = residuals(gamm_out, type = "deviance")
hake_weight_age_df$gamm_GS.resids = residuals(gamm_GS_out, type = "deviance")
hake_weight_age_df$gamm_cohort.resids = residuals(gamm_cohort_out, type = "deviance")
hake_weight_age_df$gamm_GS_cohort.resids = residuals(gamm_GS_cohort_out, type = "deviance")
hake_weight_age_df$gamm_GS_year_cohort.resids = residuals(gamm_GS_year_cohort_out, type = "deviance")
hake_weight_age_df$gamm_year.resids = residuals(gamm_year_out, type = "deviance")
hake_weight_age_df$gamm_year_cohort.resids = residuals(gamm_year_cohort_out, type = "deviance")



### brms version ------
gam_brm_out <- brm(bf(weight ~ s(new_age)),
          data = hake_weight_age_df, family = gaussian(), cores = 4,
          iter = 2000, warmup = 1000, chains = 4)
#save(gam_brm_out, file = "results/gam_brm_out.RData")
hake_weight_age_df$brm.gam.resids = residuals(gam_brm_out)[,1]


gamm_year_brm_out <- brm(bf(weight ~ s(new_age) + s(catch_year, bs = 're')),
                   data = hake_weight_age_df, family = gaussian(), cores = 4,
                   iter = 2000, warmup = 1000, chains = 4)
#save(gamm_year_brm_out, file = "results/gamm_year_brm_out.RData")
hake_weight_age_df$brm.gamm.year.resids = residuals(gamm_year_brm_out)[,1]

gamm_yr_GS_brm_out <- brm(bf(weight ~ s(new_age, m = 2) +
                               s(new_age, catch_year, bs="fs", m = 2)),
                         data = hake_weight_age_df, family = gaussian(), cores = 4,
                         iter = 2000, warmup = 1000, chains = 4)

### Save dataframe with all resids
#save(hake_weight_age_df, file = "results/hake_weight_age_df.RData")

#AIC(gamm_GS_year_cohort_out, gamm_GS_cohort_out, gamm_cohort_out,gamm_out,gamm_GS_out)

# explore resids -----------------------------------------------------


# model with different wiggliness
p2 = hake_weight_age_df %>%  #summarise residual information by year
  group_by(catch_year) %>% 
  summarise(avg = mean(gamm_cohort.resids), sd = sd(gamm_cohort.resids), n = n()) %>% 
  ggplot(aes(x = catch_year, y = avg)) +
  geom_point(aes(size = n)) +
  scale_size(range = c(1, 3)) +
  theme_classic() +
  labs(title = 'average weight-at-age anomaly', subtitle = "cohort random effect with different wiggliness", x = 'year', y = 'residual') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  lims(y = c(-0.008, 0.025))

# model with same wiggliness
p1 = hake_weight_age_df %>%  #summarise residual information by year
  group_by(catch_year) %>% 
  summarise(avg = mean(gamm_GS_cohort.resids), sd = sd(gamm_GS_cohort.resids), n = n()) %>% 
  ggplot(aes(x = catch_year, y = avg)) +
  geom_point(aes(size = n)) +
  scale_size(range = c(1, 3)) +
  theme_classic() +
  labs(title = 'average weight-at-age anomaly', subtitle = "cohort and year random effect with same wiggliness", x = 'year', y = 'residual') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  lims(y = c(-0.001, 0.001))

jpeg(file="plots/weight-at-age/avg_growth_anomaly_cohort_gams.jpeg")
ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
dev.off()

jpeg(file="plots/weight-at-age/growth_anomaly_cohort_cohortRE.jpeg")
hake_weight_age_df %>%  #summarise residual information by year
  group_by(catch_year, cohort) %>% 
  summarise(avg = mean(gamm_GS_cohort.resids), sd = sd( gamm_GS_cohort.resids), n = n()) %>% 
  filter(cohort %in% c(1996, 1999, 2005, 2010, 2013)) %>% 
  ggplot(aes(x = catch_year, y = avg, colour = as.factor(cohort), group = as.factor(cohort))) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Cohort weight-at-age anomalies through time", subtitle = "with cohort random effects", x = "Year", y = "avg growth anomaly")
dev.off()


# Condition factors -----------------------------

# cohort and year random effect
weight.age_year_cohort_gam_df = hake_weight_age_df %>% 
  select(new_age, weight, catch_year, cohort, length, distance_fished, sex_description,
         hb_latitude, hb_longitude, resids = gamm_GS_year_cohort.resids)
weight.age_year_cohort_gam_df$predictions = predict(gamm_GS_year_cohort_out)
weight.age_year_cohort_gam_df$cond.factor = weight.age_year_cohort_gam_df$weight/weight.age_year_cohort_gam_df$predictions

weight.age_year_cohort_gam_df %>%  #summarise residual information by year
  group_by(catch_year) %>% 
  summarise(avg = mean(cond.factor), sd = sd(cond.factor), n = n()) %>% 
  ggplot(aes(x = catch_year, y = avg)) +
  geom_point(aes(size = n)) +
  scale_size(range = c(1, 3)) +
  theme_classic() +
  labs(title = 'average condition factor (obs/pred)', subtitle = "cohort and year random effect with same wiggliness", x = 'year', y = 'avg condition factor') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = 1, linetype = 'dashed') 


jpeg(file="plots/weight-at-age/anomaly_spatial_peryear.jpeg")
# spatial variation
ggplot(hake_weight_age_df, aes(x = hb_longitude, y = hb_latitude, col = gamm_GS_cohort.resids)) +
  geom_point() +
  scale_colour_gradient2(limits = c(-1.42, 2.95), mid = NA) +
  theme_classic() +
  facet_wrap(vars(catch_year)) +
  scale_x_reverse()
dev.off()
# not too much visual differences between the three gam models (cohort, year, cohort + year)




# Weight at age simple EDA  ----------------

# Looking more at how residuals vary by other factors independently

# by cohorts
hake_weight_age_df %>% 
  group_by(cohort) %>% 
  summarise(avg = mean(gamm_GS_year_cohort.resids), n = n()) %>% 
  ggplot(aes(x = cohort, y = avg)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  lims(y = c(-0.05, 0.05))

# by latitude
hake_weight_age_df %>% 
  group_by(hb_latitude) %>% 
  summarise(avg = mean(gamm_GS_year_cohort.resids), n = n()) %>% 
  ggplot(aes(x = hb_latitude, y = avg)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  lims(y = c(-0.05, 0.05))

# by sex
# by latitude
hake_weight_age_df %>% 
  group_by(sex_description) %>% 
  summarise(avg = mean(gamm_GS_year_cohort.resids), n = n()) %>% 
  ggplot(aes(x = sex_description, y = avg)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  lims(y = c(-0.05, 0.05))




# plot predictions ---------
predict_data <- tidyr::expand(hake_weight_age_df, nesting(cohort, catch_year),
                          new_age = unique(new_age))

m1_pred <- bind_cols(predict_data,
                     as.data.frame(predict(gamm_year_cohort_out, newdata = predict_data,
                                           se.fit = TRUE)))

ggplot(m1_pred, aes(x = new_age, y = fit, group = catch_year, color = catch_year)) +
  geom_smooth(size = 0.1) +
  facet_wrap(vars(cohort)) +
  theme_classic()


# brms - gams ------------
# Now that I have a better understanding of gams, I will explore bayesian gams

gamm_brm_out = brm(bf(weight ~ s(new_age)), 
                        data = hake_weight_age_df, family = lognormal(), cores = 4,
                        iter = 2000, warmup = 1000, chains = 4)
#save(gamm_brm_out, file = "bayes_results/gamm_brm_out.RData")



gamm_year_brm_out = brm(bf(weight ~ s(new_age) + s(catch_year, bs="re")), 
                    data = hake_weight_age_df, family = lognormal(), cores = 4,
                    iter = 2000, warmup = 1000, chains = 4)
#msms <- conditional_smooths(gamm_year_brm_out)
#plot(msms)
save(gamm_year_brm_out, file = "bayes_results/gamm_year_brm_out.RData")

gamm_cohort_brm_out = brm(bf(weight ~ s(new_age) + s(cohort, bs="re")), data = hake_weight_age_df, 
                          family = lognormal(), cores = 4, iter = 2000, warmup = 1000, chains = 4)
#save(gamm_cohort_brm_out, file = "bayes_results/gamm_cohort_brm_out.RData")


# This is the favored one
gamm_year_cohort_out = brm(bf(weight ~ s(new_age) + s(catch_year, bs="re")
                              + s(cohort, bs="re")), data = hake_weight_age_df, family = lognormal(), 
                           cores = 4, iter = 2000, warmup = 1000, chains = 4)
#save(gamm_year_cohort_out, file = "bayes_results/gamm_year_cohort_out.RData")
msms <- conditional_smooths(gamm_cohort_brm_out)
plot(msms)

get_variables(gamm_year_cohort_out)

post_ranef <- posterior_samples(gamm_year_cohort_out, add_chain = T) %>% 
  dplyr::select(-lp__, -iter, -contains("b_"), -contains("sds_"))

mcmc_trace(post_ranef) # traceplots

rhat_vals <- rhat(gamm_year_cohort_out) #Rhat
mcmc_rhat_data(rhat_vals)
mcmc_rhat(rhat_vals) + theme_bw()

neff_vals <- neff_ratio(gamm_year_cohort_out) # NEFF
mcmc_neff_data(neff_vals)
mcmc_neff(neff_vals)  + theme_bw()

mcmc_acf(posterior_samples(gamm_year_cohort_out)) # autocorrelation


# Visualize random effect predictions
re_model_only <- tidyr::crossing(new_age = seq(min(hake_weight_age_df$new_age), 
                                        max(hake_weight_age_df$new_age), length.out=100),
                          catch_year = unique(hake_weight_age_df$catch_year),
                          cohort = unique(hake_weight_age_df$cohort)) %>%
        add_epred_draws(gamm_year_cohort_out, ndraws = 1e3)

re_model_summary <- re_model_only %>%
  group_by(catch_year, new_age) %>%
  dplyr::summarize(.epred = mean(.epred))

all_pred <- tidyr::crossing(new_age = seq(min(hake_weight_age_df$new_age), 
                                 max(hake_weight_age_df$new_age), length.out=100),
                   catch_year = unique(hake_weight_age_df$catch_year),
                   cohort = unique(hake_weight_age_df$cohort)) %>%
  add_predicted_draws(gamm_year_cohort_out,
                      allow_new_levels = TRUE,
                      ndraws = 1e3)

pp_check(gamm_year_cohort_out, ndraws = 1e2) + 
  ggtitle("PPC gamm with year+cohort RE") + 
  theme_bw(base_size = 10) +
  xlim(-1, 3)


jpeg(file="plots/weight-at-age/bayesian/predicted_weight_per_year.jpeg")
ggplot(re_model_only,
       aes(x = new_age, y = .epred)) +
  facet_wrap(~catch_year) +
  stat_interval() +
  theme_classic()
dev.off()

ggplot(re_model_only,
       aes(x = new_age, y = .epred)) +
  geom_line(aes(group = .draw), alpha = 0.1) +
  geom_line(data = re_model_summary, alpha = 0.8, 
            color = "red", lwd = 1) +
  facet_wrap(~catch_year)


# simple data exploration of weight at age trends (no modeling) -----
# following what pollock assessments did

avg_weight_at_age = hake_weight_age_df %>% 
  group_by(new_age) %>% 
  summarise(avg_weight = mean(weight))
as.data.frame(avg_weight_at_age)

hake_weight_age_df = left_join(hake_weight_age_df, avg_weight_at_age, by = "new_age")
hake_weight_age_df$std_weight_dev = hake_weight_age_df$weight/hake_weight_age_df$avg_weight
as.factor(hake_weight_age_df$new_age)

ggplot(hake_weight_age_df, aes(x = new_age, y = std_weight_dev, group = new_age)) +
  geom_violin() + 
  stat_summary(fun.data=data_summary) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = 1, lty = 2) +
  facet_wrap(vars(catch_year))



# To compare models

loo(gamm_year_cohort_out, gamm_year_brm_out, gamm_cohort_brm_out)



# let's figure out what the peaks are
ggplot(hake_weight_age_df, aes(x = weight, group = as.factor(cohort), fill = as.factor(cohort))) +
  geom_histogram(bins = 60) +
  theme_classic() +
  scale_fill_discrete(limits = c(0,1,2,3,4,5,6,7))


# testing different distributions using mgcv

gamm_cohort_gamma_out = gam(weight ~ s(new_age) + s(cohort), data = hake_weight_age_df, 
                            family = Gamma())
gamm_cohort_lnorm_out = gam(weight ~ s(new_age) + s(cohort, bs="re"), data = hake_weight_age_df, 
                            family = gaussian(link = "log"))


gamm_cohort_year_gamma_out = gam(weight ~ s(new_age, k = 12) + s(catch_year, bs="re") + s(cohort, k = 20), data = hake_weight_age_df, 
                            family = Gamma())
gamm_cohort_year_lnorm_out = gam(weight ~ s(new_age, k = 12) + s(catch_year, bs="re") + s(cohort, k = 17), data = hake_weight_age_df, 
                            family = gaussian(link = "log"))

gamm_cohort_year_space_lnorm_out = gam(weight ~ s(new_age, k = 10) + s(catch_year, bs="re") + s(cohort, k = 18) +
                                   s(hb_latitude, k = 20), data = hake_weight_age_df, 
                                 family = gaussian(link = "log"))

rsd = residuals(gamm_cohort_year_space_lnorm_out)
prd = predict(gamm_cohort_year_space_lnorm_out)
hake_weight_age_df = hake_weight_age_df %>% 
  mutate(resids = rsd, preds = exp(prd))

ggplot(hake_weight_age_df, aes(x = weight, y = resids)) +
  geom_point() +
  facet_wrap(vars(new_age))



# Let's check correlations....
hake_corr_check = hake_weight_age_df %>% 
  select(catch_month, catch_day, length, distance_fished, hb_latitude, hb_longitude, catch_year, cohort)
hake_corr_check$cohort = as.numeric(hake_corr_check$cohort)
hake_corr_check$catch_year = as.numeric(hake_corr_check$catch_year)
cor(hake_corr_check)


# plotting weight at age density plot

ggplot(hake_weight_age_df, aes(x=weight)) + 
  geom_density() +
  lims(x = c(0,1.5)) +
  theme_classic() +
  facet_wrap(~catch_year, ncol = 2,  dir = "v") +                                          
  geom_vline(data = median_weights_by_year, aes(xintercept = med_weight), linetype = "dashed", col = "blue")

median_weights_by_year = as.data.frame(hake_weight_age_df %>% 
  group_by(catch_year) %>% 
  summarise(med_weight = median(weight), n = n()))


