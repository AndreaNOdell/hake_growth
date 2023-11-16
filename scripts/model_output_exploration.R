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
library(stringr) 


# Survey --------
load("results/hake_weight_age_df_updated.RData")
load("results/sdmTMB/m_age_month_cohort_sex_insmooth.RData")
model = m_age_month_cohort_sex_insmooth
df = model$data

# model diagnostics
coefs = rbind(tidy(model, effects = "ran_pars", conf.int = TRUE), tidy(model, effects = "fixed", conf.int = TRUE))
df$residuals <- residuals(model) # randomized quantile residuals
predict_df = predict(model)
df = cbind(df, predict_df[,(ncol(predict_df)-3):ncol(predict_df)])


## Residuals ---------------
jpeg(filename = "plots/output_exploration/residuals_peryear_spatial.jpeg", units="in", width=5, height=5, res = 300)
year_sampled = sort(unique(hake_weight_age_df_updated$catch_year))
ggplot(df[df$catch_year %in% year_sampled,], aes(X, Y, col = residuals)) +
  scale_colour_gradient2() +
  geom_point() +
  #facet_wrap(~catch_year) +
  coord_fixed() +
  labs(title= "model residuals") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme_classic()
dev.off()

## Data trends in weight by age --------------
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

## Data trends by sex -----------
hake_sdmTMB_df %>% 
  group_by(sex_description) %>% 
  summarise(n = n()) # 138 observations unsexed out of 35,348 (< 1%) and about 90% of these observations are in early age classes (0 and 1)

## Point estimates ---------
jpeg(file="plots/updated_output_exploration/point_estimates.jpeg", units="in", width=6.5, height=4, res = 300)
coefs %>% 
  filter(!term %in% "range") %>% 
  ggplot(aes( x = estimate, y = term)) +
  geom_point() +
  theme_classic() +
  lims( x = c(-1.65, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high), width=.1) +
  labs(y = "coefficient", subtitle = "point estimates with confidence intervals")
dev.off()

options(knitr.kable.NA = '')
kbl(coefs, align = "c") %>%
  kable_classic_2(full_width = F) %>% 
  #add_header_above(data.frame("centered data",6), monospace = TRUE) %>%
save_kable(file = "plots/updated_output_exploration/coefficient_table.jpeg", zoom = 4)


## Cohort Random Effects intercepts ----------
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

## cohort smoother ----------
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


## Predicted weight by year -----------
jpeg(filename = "plots/updated_output_exploration/predicted_weight_year.jpeg", units="in", width=4, height=3, res = 300)
pred_out %>% 
  filter(catch_year %in% year_sampled) %>% 
  group_by(new_age, catch_year) %>% 
  summarise(avg = mean(exp(est)), sd = sd(exp(est))) %>% 
  filter(new_age %in% 2:10) %>% 
  ggplot(aes(x = catch_year, y = avg, group = as.factor(new_age))) +
    labs(y = "Predicted weight-at-age", x = "Year") +
    #geom_ribbon(aes(ymin = avg - sd, ymax = avg + sd, fill = as.factor(new_age)), alpha = 0.5) +
    geom_line(aes(col = as.factor(new_age))) +
  theme_classic()
dev.off()

year_sample_size = hake_weight_age_df_updated %>% 
  group_by(catch_year, new_age) %>% 
  summarise(n = n())

ggplot(year_sample_size, aes(x = catch_year, y = n, group = as.factor(new_age), fill = as.factor(new_age))) +
  geom_col() +
  theme_classic() 

# cohort effects
pred_out %>% 
  filter(new_age %in% c(5, 10)) %>% 
  group_by(cohort, new_age) %>% 
  summarise(avg = mean(exp(est)), sd = sd(exp(est))) %>% 
  ggplot(aes(x = cohort, y = avg, group = as.factor(new_age))) +
  labs(y = "Predicted weight at age", x = "Cohort") +
  geom_ribbon(aes(ymin = avg - sd, ymax = avg + sd, fill = as.factor(new_age)), alpha = 0.5) +
  geom_line(aes(col = as.factor(new_age))) +
  theme_classic() +
  geom_vline(xintercept = c(1980, 1984, 1999, 2010), linetype = "dashed", col = "gray") 

pred_out %>% 
  group_by(cohort, new_age) %>% 
  summarise(avg = mean(exp(est)), sd = sd(exp(est))) %>% 
  filter(cohort %in% c(1984:2015)) %>% 
  ggplot(aes(x = new_age, y = avg, group = cohort)) +
  labs(y = "Predicted weight at age", x = "Age") +
  #geom_ribbon(aes(ymin = avg - sd, ymax = avg + sd, fill = as.factor(cohort)), alpha = 0.5) +
  geom_line(aes(col = cohort)) +
  theme_classic()
  

pred_out %>% 
  filter(catch_year %in% c(2000, 2003, 2009, 2011, 2013, 2015, 2017)) %>% 
  group_by(cohort, catch_year) %>% 
  summarise(avg = mean(exp(est)), sd = sd(exp(est))) %>% 
  #filter(new_age %in% c(3, 5, 10)) %>% 
  ggplot(aes(x = cohort, y = avg, group = as.factor(catch_year))) +
  labs(y = "Predicted weight at age", x = "Cohort") +
  geom_ribbon(aes(ymin = avg - sd, ymax = avg + sd, fill = as.factor(catch_year)), alpha = 0.5) +
  geom_line(aes(col = as.factor(catch_year))) +
  theme_classic() +
  geom_vline(xintercept = c(1980, 1984, 1999, 2010))


# Year - age = cohort
# Year - cohort = Age

# Plot 2D effects of age and year -- as anomalies
pred_out <- group_by(pred_out, sex_description, new_age) %>% 
  dplyr::mutate(m = mean(est),
                anom_est = est-m)

anom <- pred_out %>%
  mutate(rem = catch_year-cohort) %>%
  filter(rem > 0) %>%
  dplyr::rename(Anomaly = anom_est) %>% 
  group_by(cohort, catch_year) %>% 
  summarise(avg_anom = mean(Anomaly), n = n())

anom_fem <- dplyr::filter(pred_out,sex_description=="Female") %>%
  mutate(rem = catch_year-cohort) %>%
  filter(rem > 0) %>%
  dplyr::rename(Anomaly = anom_est) %>% 
  group_by(cohort, catch_year) %>% 
  summarise(avg_anom = mean(Anomaly), n = n())

anom_male <- dplyr::filter(pred_out,sex_description=="Male") %>%
  mutate(rem = catch_year-cohort) %>%
  filter(rem > 0) %>%
  dplyr::rename(Anomaly = anom_est) %>% 
  group_by(cohort, catch_year) %>% 
  summarise(avg_anom = mean(Anomaly), n = n())

jpeg(filename = "plots/updated_output_exploration/Year_Cohort_Anomalies.jpeg", units="in", width=8, height=5, res = 300)
ggplot(anom, aes(catch_year, cohort, col=avg_anom)) + 
  geom_point(aes(size = n)) + 
  labs(x = "Year effect", y = "Cohort effect", col = "Anomaly", size = "sample size") + 
  theme_classic() + 
  scale_color_gradient2(limits = c(-0.35,0.47)) + 
  theme(strip.background=element_rect(fill="white")) 
dev.off()


# do heavier individuals in specific age-classes occur more often in the north?

pred_out %>% 
  filter(new_age %in% 5) %>% 
  group_by(catch_month, new_age) %>% 
  summarise(avg = mean(exp(est)), sd = sd(exp(est))) %>% 
  ggplot(aes(x = catch_month, y = avg)) +
  labs(y = "Predicted weight of age 5 individual", x = "month") +
  lims(y = c(0.35,0.8)) +
  geom_line() +
  theme_classic()

pred_out %>% 
  group_by(sex_description, new_age) %>% 
  summarise(avg = mean(exp(est)), n = n()) %>% 
  ggplot(aes(x = new_age, y = avg, group = sex_description)) +
  geom_point(aes(size = n, col = sex_description)) +
  geom_line(aes(col = sex_description)) +
  labs(x = "age", y = "avg weight") +
  theme_classic()


# Fishery ----------
load("results/sdmTMB/fm_age_sex_cohort_month_yearre.RData") # load in data
model_fishery = fm_age_sex_cohort_month_yearre   # assign model a name
df_fishery = model_fishery$data # assign data a name

coefs_fishery = rbind(tidy(model_fishery, effects = "ran_pars", conf.int = TRUE), tidy(model_fishery, effects = "fixed", conf.int = TRUE), tidy(model_fishery, effects = "ranef", conf.int = TRUE))
df_fishery$residuals <- residuals(model_fishery) # randomized quantile residuals
predict_df_fishery = predict(model_fishery)
df_fishery = cbind(df_fishery, est = predict_df_fishery[,ncol(predict_df_fishery)])


View(df_fishery %>% 
  group_by(catch_month, catch_year,country) %>% 
  summarise(n = n()) %>% 
  filter(catch_month %in% c(12, 1, 2, 3)))

df_fishery %>% 
  filter(new_age %in% c(10)) %>% 
  group_by(country, new_age, catch_year) %>% 
  summarise(avg = mean(exp(est)), n = n()) %>% 
  ggplot(aes(x = catch_year, y = avg, group = country)) +
  geom_point(aes(size = n, col = country)) +
  geom_line(aes(col = country)) +
  labs(x = "year", y = "avg weight of 10 y/o") +
  theme_classic()

df_fishery %>% 
  group_by(sex_description, new_age) %>% 
  summarise(avg = mean(exp(est)), n = n()) %>% 
  ggplot(aes(x = new_age, y = avg, group = sex_description)) +
  geom_point(aes(size = n, col = sex_description)) +
  geom_line(aes(col = sex_description)) +
  labs(x = "age", y = "avg weight") +
  theme_classic()

df_fishery %>% 
  group_by(country, new_age) %>% 
  summarise(avg = mean(exp(est)), n = n()) %>% 
  ggplot(aes(x = new_age, y = avg, group = country)) +
  geom_point(aes(size = n, col = country)) +
  geom_line(aes(col = country)) +
  labs(x = "age", y = "avg weight") +
  theme_classic()

# Combined dataframes ----------
fishery_for_combined = df_fishery %>% 
  select("sex_description", "catch_month", "catch_year", "new_age", "cohort", "country", "fcatch_year", "est") %>% 
  mutate(data_source = "fishery") %>% 
  mutate(across('sex_description', str_replace, 'U', 'Undetermined')) %>% 
  mutate(across('sex_description', str_replace, 'F', 'Female')) %>% 
  mutate(across('sex_description', str_replace, 'M', 'Male')) 
  
fishery_for_combined$new_age = as.integer(fishery_for_combined$new_age)

survey_for_combined = as.data.frame(pred_out) %>% 
  select("sex_description", "catch_month", "catch_year", "new_age", "cohort", "est", "hb_latitude") %>% 
  mutate(fcatch_year = as.factor(catch_year), country = as.factor(ifelse(hb_latitude > 49, "canada", "usa")), data_source = "survey") %>% 
  select("sex_description", "catch_month", "catch_year", "new_age", "cohort", "country", "fcatch_year", "est", "data_source")
survey_for_combined$catch_month = as.integer(survey_for_combined$catch_month)
survey_for_combined$catch_year = as.integer(survey_for_combined$catch_year)

combined_df = rbind(survey_for_combined, fishery_for_combined)
combined_df$data_source = as.factor(combined_df$data_source)
#save(combined_df, file = "results/sdmTMB/combined_df.RData")


jpeg(filename = "plots/updated_output_exploration/combined/10yoWeightTimeseries_DataSource.jpeg", units="in", width=8, height=5, res = 300) 
combined_df %>% 
  filter(new_age %in% c(10)) %>% 
  group_by(data_source, new_age, catch_year) %>% 
  summarise(avg = mean(exp(est)), n = n()) %>% 
  ggplot(aes(x = catch_year, y = avg, group = data_source)) +
  geom_point(aes(size = n, col = data_source)) +
  geom_line(aes(col = data_source)) +
  labs(x = "year", y = "avg weight of 10 y/o") +
  theme_classic()
dev.off()

combined_df %>% 
  filter(new_age %in% c(10)) %>% 
  group_by(country, new_age, catch_year) %>% 
  summarise(avg = mean(exp(est)), n = n()) %>% 
  ggplot(aes(x = catch_year, y = avg, group = country)) +
  geom_point(aes(col = country)) +
  geom_line(aes(col = country)) +
  labs(x = "year", y = "avg weight of 10 y/o") +
  theme_classic()

combined_df %>% 
  filter(new_age %in% c(5)) %>% 
  group_by(data_source, new_age, catch_year) %>% 
  summarise(avg = mean(exp(est)), n = n()) %>% 
  ggplot(aes(x = catch_year, y = avg, group = data_source)) +
  geom_point(aes(size = n, col = data_source)) +
  geom_line(aes(col = data_source)) +
  labs(x = "year", y = "avg weight of 5 y/o") +
  theme_classic()
# the survey and fishery trends are fairly similar, with the exception of a few years (1993-1996 and 2007ish). We can see from the size of the points the sample size from each data source for the year. The survey contributes are larger sample size, but only for a subset of the years.

jpeg(filename = "plots/updated_output_exploration/combined/weight-at-age_sex_MFonly.jpeg", units="in", width=12, height=7, res = 600) 
combined_df %>% 
  filter(catch_year >= min(survey_for_combined$catch_year)) %>%  
  filter(!sex_description == "Undetermined") %>% 
  group_by(sex_description, new_age) %>%  
  summarise(avg = mean(exp(est)), n = n(), sd = sd(exp(est))) %>% 
  ggplot(aes(x = new_age, y = avg, ymin = avg-sd, ymax = avg+sd, linetype = sex_description)) +
  #geom_point(aes(shape = sex_description)) + 
  geom_line(size = 1.3) + 
  geom_ribbon(alpha = 0.25) +
  scale_color_manual(values=c("#F05039", "#1F449C")) +
  labs(x = "age", y = "estimated weight", linetype = "sex") +
  theme_classic() +
  lims(y = c(0,1.2)) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))
dev.off()

jpeg(filename = "plots/updated_output_exploration/combined/weight-at-age_sex_data-source_MFonly.jpeg", units="in", width=12, height=7, res = 600)
combined_df %>% 
  filter(catch_year >= min(survey_for_combined$catch_year)) %>% 
  filter(!sex_description == "Undetermined") %>% 
  group_by(data_source, sex_description, new_age) %>%  
  summarise(avg = mean(exp(est)), n = n(), sd = sd(exp(est))) %>% 
  ggplot(aes(x = new_age, y = avg)) +
  #geom_point(aes(col = data_source, shape = sex_description)) + 
  geom_line(aes(col = data_source, linetype = sex_description), size = 1.3) + 
  scale_color_manual(values=c("#F05039", "#1F449C")) +
  labs(x = "age", y = "estimated weight", linetype = "sex", color = "data source") +
  theme_classic()+
    lims(y = c(0,1)) +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18)) +
  guides(linetype = "none")
# the weight representation is similar between the fishery and survey for both males and females. The undetermined sex, which is primarily the fishery in Canada, has greater weight at age than the U.S.
dev.off()


jpeg(filename = "plots/updated_output_exploration/combined/weight-at-age_data-source.jpeg", units="in", width=8, height=5, res = 300)
combined_df %>% 
  filter(catch_year >= min(survey_for_combined$catch_year)) %>%  
  group_by(data_source, new_age) %>% 
  summarise(avg = mean(exp(est)), n = n()) %>% 
  ggplot(aes(x = new_age, y = avg)) +
  geom_point(aes(col = data_source, shape = data_source)) +
  geom_line(aes(col = data_source, linetype = data_source)) +
  labs(x = "age", y = "avg weight") +
  theme_classic()
dev.off()

jpeg(filename = "plots/updated_output_exploration/combined/weight-at-age_sex_country.jpeg", units="in", width=12, height=7, res = 600)
combined_df %>% 
  filter(catch_year >= min(survey_for_combined$catch_year)) %>%  
  group_by(sex_description, country, new_age) %>% 
  summarise(avg = mean(exp(est)), n = n()) %>% 
  ggplot(aes(x = new_age, y = avg)) +
  #geom_point(aes(col = sex_description, shape = country)) +
  geom_line(aes(col = sex_description, linetype = country), size = 1) +
  scale_color_manual(values=c("#F05039", "#1F449C", "darkgray")) +
  labs(x = "age", y = "estimated weight", col = "sex") +
  theme_classic() +
  lims(y = c(0,1)) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)) 
dev.off()

jpeg(filename = "plots/updated_output_exploration/combined/weight-at-age_country_datasource.jpeg", units="in", width=12, height=7, res = 600)
combined_df %>% 
  filter(catch_year >= min(survey_for_combined$catch_year)) %>% 
  #filter(!sex_description == "Undetermined") %>% 
  group_by(data_source, country, new_age) %>% 
  summarise(avg = mean(exp(est)), n = n()) %>% 
  ggplot(aes(x = new_age, y = avg)) +
  #geom_point(aes(col = sex_description, shape = country)) +
  geom_line(aes(linetype = country, color = data_source), size = 1) + 
  scale_color_manual(values=c("#F05039", "#1F449C")) +
  labs(x = "age", y = "estimated weight", col = "data source") +
  theme_classic() +
  lims(y = c(0,1)) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18)) +
  guides(linetype = "none")
dev.off()

jpeg(filename = "plots/updated_output_exploration/combined/weight-at-age_country.jpeg", units="in", width=12, height=7, res = 600)
combined_df %>% 
  filter(catch_year >= min(survey_for_combined$catch_year)) %>% 
  #filter(!sex_description == "Undetermined") %>% 
  group_by(country, new_age) %>% 
  summarise(avg = mean(exp(est)), n = n()) %>% 
  ggplot(aes(x = new_age, y = avg)) +
  #geom_point(aes(col = sex_description, shape = country)) +
  geom_line(aes(linetype = country), size = 1) + 
  scale_color_manual(values=c("#6e7cb9", "#d2848d")) +
  labs(x = "age", y = "estimated weight", col = "data source") +
  theme_classic() +
  lims(y = c(0,1)) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))
dev.off()


# anomalies 
combined_df <- group_by(combined_df, data_source, sex_description, new_age) %>% 
  dplyr::mutate(m = mean(est),
                anom_est = est-m)

jpeg(filename = "plots/updated_output_exploration/combined/year_cohort_anomalies_wide.jpeg", units="in", width=14, height=7, res = 300)
combined_df %>%
  #filter(catch_year >= min(survey_for_combined$catch_year)) %>%  
  mutate(rem = catch_year-cohort) %>%
  filter(rem > 0) %>%
  dplyr::rename(Anomaly = anom_est) %>% 
  group_by(data_source, cohort, catch_year) %>% 
  summarise(avg_anom = mean(Anomaly), n = n()) %>% 
  #filter(cohort %in% c(1982, 1999)) %>% 
ggplot(aes(catch_year, cohort, col=avg_anom)) + 
  geom_point(size = 6) + 
  labs(x = "Year effect", y = "Cohort effect", col = "anomaly") + #, size = "sample size" 
  theme_classic() + 
  scale_color_gradient2(high = "green4", mid = "white", low = "darkorchid3",limits = c(-0.626, 0.915)) + 
  theme(strip.background=element_rect(fill="white"),
        #strip.text.x = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 20))+
  lims(y = c(min(combined_df$cohort), max(combined_df$cohort)),
       x = c(min(combined_df$catch_year), max(combined_df$catch_year))) +
  facet_wrap(~data_source, ncol = 1) #+
  #geom_hline(yintercept = c(1980, 1984, 1999, 2010, 2014), linetype = 2, col = "gray")
dev.off()


# I collected recruitment proportions for pacific hake from table 10 of the 2021 stock assessment to see whether there was a correlation between large recruitment events and impacts on weight-at-age (i.e. cohort effects)
recruit_props = as.data.frame(cbind(recr_prop = c(4.61, 0.09, 0.00, 0.47, 0.00, 0.15, 19.49, 0,0,0,0.92,0,0,0,0,0,0,0.46,0,0,0.62,0,0,0.02,0.06,1,0,0,0,0,0.02,0.33,0.78,0.76,0.64,0.03,2.67,0.18,0.03,0,3.64,0.29,3.76,7.35,0.01,0), cohort = 1974:2019))

combined_df = as.data.frame(left_join(combined_df, recruit_props))


combined_df %>% 
  filter(new_age == 4) %>%
  filter(!is.na(recr_prop)) %>% 
  #filter(recr_prop >= 1) %>% 
  ggplot(aes(x = recr_prop, y = anom_est)) +
  geom_point() +
  labs(y = "Predicted weight anomaly of 12y/o", x = "Recruitment Intensity") +
  theme_classic() +
  geom_hline(yintercept = 0, col = "gray", linetype = 2)


x = combined_df %>% 
  filter(new_age %in% c(3,12)) %>% 
  group_by(cohort, new_age) %>% 
  summarise(avg = mean(anom_est))
  
rem = x$cohort[!(duplicated(x$cohort)|duplicated(x$cohort,fromLast=TRUE))]
x = x[!x$cohort %in% rem, ]

x = x %>% 
  pivot_wider(names_from = new_age, values_from = avg) %>% 
  mutate(diff = (`12` - `3`)/`3`)
x = left_join(x, recruit_props) %>% 
  filter(!is.na(recr_prop)) %>% 
  mutate(year = cohort + 12)

ggplot(x, aes(x = recr_prop, y = diff)) +
  geom_point() +
  theme_classic() +
  geom_hline(yintercept = 0, col = "gray", linetype = 2) +
  labs(x = "recruitment intensity", y = "change in avg weight anomaly \n between 3 and 12 y/o") # this plot indicates that cohort effects are not necessarily driven by recruitment magnitude 

ggplot(x, aes(x = year, y = diff)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  geom_hline(yintercept = 0, col = "gray", linetype = 2) +
  labs(x = "year", y = "change in avg weight anomaly \n between 3 and 12 y/o")

x_all = combined_df %>% 
  group_by(cohort, new_age) %>% 
  summarise(avg = mean(anom_est))
rem_xall = x_all$cohort[!(duplicated(x_all$cohort)|duplicated(x_all$cohort,fromLast=TRUE))]
x_all = x_all[!x_all$cohort %in% rem, ]
x_all = x_all %>% 
  pivot_wider(names_from = new_age, values_from = avg)

# sample size across years and data sources
combined_df %>% 
  group_by(catch_year, data_source) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = catch_year, y = n, fill = data_source)) +
  geom_col(aes(position = "dodge"))

# sample size by age/sex for fishery data
jpeg(filename = "plots/updated_output_exploration/combined/fishery_samplesize_age.jpeg", units="in", width=8, height=6, res = 300)
combined_df %>% 
  filter(data_source == "fishery") %>% 
  group_by(catch_year, new_age) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = catch_year, y = n)) +
  geom_col(aes(fill = new_age)) 
dev.off()


jpeg(filename = "plots/updated_output_exploration/combined/cohort_anomalies_timeseries.jpeg", units="in", width=8, height=5, res = 300)
combined_df %>% 
  #filter(cohort %in% c(1996, 1999, 2004, 2007, 2009)) %>% 
  group_by(catch_year, cohort) %>% 
  summarise(avg = mean(anom_est)) %>% 
  ggplot(aes(x = catch_year, y = avg, group = cohort)) +
  geom_point(aes(col = cohort), size = 2) +
  geom_line(aes(col = cohort), size = 0.75) +
  theme_classic() +
  geom_hline(yintercept = 0, col = "gray", linetype = 2) +
  labs(x = "Year", y = "Predicted weight-at-age anomaly", col = "Cohort") +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14))
dev.off()

# Year effects seem to outweight cohort effects.


df_fishery %>% 
  filter(new_age == 10) %>% 
  ggplot(aes(x = catch_month, y = weight)) +
  geom_point()


