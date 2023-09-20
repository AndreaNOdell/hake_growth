library(dplyr)
library(sf)
library(sp)
library(rgdal)
library(ggplot2)
library(sdmTMB)
library(visreg)
library(ggpubr)
library(ggbreak)
library(kableExtra)

load("results/hake_weight_age_df.RData")

load("results/sdmTMB/m_age_month_cohort.RData")
load("results/sdmTMB/cohortsmoother_s_st_model.RData")


model = m_age_month_cohort
df = model$data

# model diagnostics
coefs = rbind(tidy(model, effects = "ran_pars", conf.int = TRUE), tidy(model, effects = "fixed", conf.int = TRUE))
df$residuals <- residuals(model) # randomized quantile residuals
predict_df = predict(model)
df = cbind(df, predict_df[,(ncol(predict_df)-3):ncol(predict_df)])


# Residuals ---------------
jpeg(filename = "plots/output_exploration/residuals_peryear_spatial.jpeg", units="in", width=5, height=5, res = 300)
ggplot(df[df$catch_year %in% year_sampled,], aes(X, Y, col = residuals)) +
  scale_colour_gradient2() +
  geom_point() +
  #facet_wrap(~catch_year) +
  coord_fixed() +
  labs(title= "model residuals") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme_classic()
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

# Data trends by sex -----------
hake_sdmTMB_df %>% 
  group_by(sex_description) %>% 
  summarise(n = n()) # 138 observations unsexed out of 35,348 (< 1%) and about 90% of these observations are in early age classes (0 and 1)

# Point estimates ---------
jpeg(file="plots/nested_models/cohort_month_sex_model/point_estimates.jpeg", units="in", width=6.5, height=4, res = 300)
coefs %>% 
  filter(!term %in% "range") %>% 
  ggplot(aes( x = estimate, y = term)) +
  geom_point() +
  theme_classic() +
  lims( x = c(-2, 0.82)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high), width=.1) +
  labs(y = "coefficient", subtitle = "point estimates with confidence intervals")
dev.off()

options(knitr.kable.NA = '')
kbl(coefs, align = "c") %>%
  kable_classic_2(full_width = F) %>% 
  #add_header_above(data.frame("centered data",6), monospace = TRUE) %>%
save_kable(file = "plots/output_exploration/coefficient_table.jpeg", zoom = 4)


# Cohort Random Effects intercepts ----------
year_sampled = as.numeric(sort(unique(hake_weight_age_df$catch_year))) # specify extra time vector
varying_intercepts = tidy(m_age_cohort_re_s_st, "ranef", conf.int = TRUE)[,c(1:2,4:5)] # extract random intercepts
varying_intercepts$term <- c(1973:2016) # rename years
varying_intercepts = varying_intercepts %>% 
  mutate(year_sampled = term %in% year_sampled)

jpeg(filename = "plots/nested_models/cohortRE_model/cohort_random_effects_s_st.jpeg", units="in", width=7, height=3, res = 300)
ggplot(varying_intercepts, aes(x = term, y = estimate)) +
  geom_errorbar(inherit.aes = FALSE, aes(x = term, ymax = conf.high, ymin = conf.low, 
                                        color = year_sampled), size = 0.2) +
  geom_point(aes(color = year_sampled)) +
  scale_color_manual(values = c("gray", "black")) +
  theme_classic() +
  labs(title = "cohort s-st model random effects")
dev.off()

# cohort smoother ----------
jpeg(filename = "plots/output_exploration/cohort_GAM.jpeg", units="in", width=5, height=4, res = 300)
visreg::visreg(model, xvar = "cohort", xlim = c(1980, 2017), ylim = c(-5, 5), data = model$data)
dev.off()

jpeg(filename = "plots/nested_models/cohort_smoother_model/cohort-GAM-scaled.jpeg", units="in", width=4, height=3, res = 300)
visreg::visreg(m_age_cohort_sm_s_st, xvar = "cohort", 
               xlim = c(1980, 2015), ylim = c(-1, 2), scale = "response")
#title(main = "cohort smoothed function")
dev.off()

jpeg(filename = "plots/output_exploration/age_GAM.jpeg", units="in", width=5, height=4, res = 300)
visreg::visreg(model, xvar = "new_age", xlim = c(0, 15), ylim = c(-7, 7))
dev.off()

jpeg(filename = "plots/nested_models/cohort_smoother_model/age-GAM-scaled.jpeg", units="in", width=4, height=3, res = 300)
visreg::visreg(model, xvar = "new_age", 
               xlim = c(0, 15), ylim = c(-0.5, 3), scale = "response")
dev.off()


# Predicted weight by year -----------
predict_avg_sim = predict(model, return_tmb_object = TRUE)
index_avg_sim <- get_index(predict_avg_sim, bias_correct = TRUE)
ggplot(index_avg_sim, aes(x = catch_year, y = est)) +
  labs(y = "Predicted weight", x = "Year") +
  geom_point(size = 2.5, shape = 21, fill = "gray1", color = "white") +
  geom_errorbar(data = index_avg_sim, inherit.aes = FALSE, 
                aes(x = catch_year, ymax = upr, ymin = lwr), alpha = 0) +
  theme(axis.title = element_text(size = 9),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.title = element_text(size = 9),
        legend.spacing.x = unit(0, 'cm'),
        plot.margin = unit(c(0.4, 0.4, 0.4, 0), "cm")) +
  theme_classic() +
  lims(y = c(0,2000))
  #guides(linetype = "none",
  #       color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1,
  #                            keywidth = 0.1, keyheight = 0.1, default.unit = "inch",
  #                            override.aes = list(linetype = c(1,1,1,1,1,0),
  #                                                shape = c(NA,NA,NA,NA,NA,16),
  #                                                alpha = rep(0.8, 6),
  #                                                color = pal)))




est <- as.list(model$sd_report, "Est")
