library(tidyverse)
library(brms)
library(mgcv)
library(ggpubr)

df = read.csv("raw_data/selection.csv") # load in data
#colnames(df)

hake_df = df %>% # make new dataframe with relevant information only
  filter(scientific_name == "Merluccius productus") %>% 
  select("age", "avg_weight", "distance_fished", #"hb_date", 
         "length", "sex_description", "weight", "catch_modified_date",
         "eq_date", "hb_date", 'hb_latitude', 'hb_longitude') %>% 
  separate(col = eq_date , into = c('date', 'time'), sep = " ", remove = FALSE) %>% 
  separate(col = date, into = c("catch_year", "catch_month", 'catch_day'), sep = "-")

# Set missing catch_year values to the year inputed for the other date related columns
hake_df[is.na(hake_df$catch_month),]$catch_year <- 2017

hake_df = hake_df %>% 
  mutate(new_age = age)
hake_df$new_age = as.integer(lapply(hake_df$new_age, function(x) ifelse(x > 14, 15, x)))
hake_df$cohort = as.numeric(hake_df$catch_year) - hake_df$new_age

ggplot(hake_df, aes(x = new_age, y = weight)) +
  geom_point() +
  theme_classic()


# Create a dataframe with age and weight, complete cases only
hake_weight_age_df = hake_df %>% 
  select(new_age, weight, catch_year, cohort, length, distance_fished, sex_description, hb_latitude, hb_longitude) %>% 
  drop_na() 

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
  labs(title = 'average weight-at-age anomaly', subtitle = "cohort random effect with same wiggliness", x = 'year', y = 'residual') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  lims(y = c(-0.008, 0.025))

jpeg(file="plots/weight-at-age/avg_growth_anomaly_cohort_gams.jpeg")
ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
dev.off()

jpeg(file="plots/weight-at-age/growth_anomaly_cohort_cohortRE.jpeg")
hake_weight_age_df %>%  #summarise residual information by year
  group_by(catch_year, cohort) %>% 
  summarise(avg = mean(gamm_cohort.resids), sd = sd( gamm_cohort.resids), n = n()) %>% 
  filter(cohort %in% c(1996, 1999, 2005, 2010, 2013)) %>% 
  ggplot(aes(x = catch_year, y = avg, colour = as.factor(cohort), group = as.factor(cohort))) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Cohort weight-at-age anomalies through time", subtitle = "with cohort random effects", x = "Year", y = "avg growth anomaly")
dev.off()

