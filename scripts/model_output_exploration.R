library(dplyr)
library(sf)
library(sp)
library(rgdal)
library(ggplot2)
library(sdmTMB)
library(visreg)
library(ggpubr)
load("results/sdmTMB/cohortRE_s_st_model.RData")
df = m_age_cohort_re_s_st$data

load("results/sdmTMB/cohortsmoother_s_st_model.RData")
df = m_age_cohort_sm_s_st$data

# model diagnostics
coefs = tidy(m_age_cohort_re_s_st, "ran_pars", conf.int = TRUE)
df$residuals <- residuals(m_age_cohort_re_s_st) # randomized quantile residuals
predict_df = predict(m_age_cohort_re_s_st)
df = cbind(df, predict_df[,(ncol(predict_df)-4):ncol(predict_df)])

# Residuals ---------------
jpeg(filename = "plots/nested_models/residuals_spatial.jpeg", units="in", width=5, height=5, res = 300)
ggplot(df, aes(X, Y, col = residuals)) +
  scale_colour_gradient2() +
  geom_point() +
  #facet_wrap(~catch_year) +
  coord_fixed() +
  labs(title= "model residuals") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# Data trends in weight by age --------------
# whole spatial domain
jpeg(filename = "plots/data_exploration/weight_trend_by_age.jpeg", units="in", width=6, height=4, res = 300)
df %>% 
  group_by(new_age, catch_year) %>% 
  summarise(avg_weight = mean(weight)) %>% 
  ggplot(aes(x = catch_year, y = avg_weight, group = as.factor(new_age), col = as.factor(new_age))) +
  geom_line() +
  labs(title = "weight through time by age", y = "avg weight", x = "year") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# differences between Canada and US
north_df = subset(df, hb_latitude > 49)
south_df = subset(df, hb_latitude < 49)
p1 = north_df %>% 
  group_by(new_age, catch_year) %>% 
  summarise(avg_weight = mean(weight)) %>% 
  ggplot(aes(x = catch_year, y = avg_weight, group = as.factor(new_age), col = as.factor(new_age))) +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Canada", y = "avg weight", x = "year") +
  lims(y = c(0,4))
p2 = south_df %>% 
  group_by(new_age, catch_year) %>% 
  summarise(avg_weight = mean(weight)) %>% 
  ggplot(aes(x = catch_year, y = avg_weight, group = as.factor(new_age), col = as.factor(new_age))) +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "USA", y = "avg weight", x = "year") +
  lims(y = c(0,4))

jpeg(filename = "plots/data_exploration/weight_trend_by_AgexCountry.jpeg", units="in", width=6, height=4, res = 300)
ggarrange(p1,p2,labels = c("A", "B"), common.legend = TRUE, legend = "right")
dev.off()


jpeg(file="plots/data_exploration/observations_of_age_by_month.jpeg", units="in", width=5, height=4, res = 300)
df %>% 
  group_by(catch_month, new_age) %>% 
  summarise(avg_weight = mean(weight), n = n()) %>% 
  ggplot(aes(x = catch_month, y = n, group = as.factor(new_age), fill = as.factor(new_age))) +
  geom_col() +
  theme_classic() 
dev.off()

# Point estimates ---------
jpeg(file="plots/nested_models/point_estimates.jpeg", units="in", width=6.5, height=4, res = 300)
coefs %>% 
  filter(!term %in% "range") %>% 
  ggplot(aes( x = estimate, y = term)) +
  geom_point() +
  theme_classic() +
  lims( x = c(-0.05, 0.62)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high), width=.1) +
  labs(y = "coefficient", subtitle = "point estimates with confidence intervals")
dev.off()


