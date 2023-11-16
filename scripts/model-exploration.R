load("results/hake_weight_age_df_updated.RData")
library(dplyr)
library(sf)
library(sp)
library(rgdal)
library(ggplot2)
library(sdmTMB)
library(visreg)

# Fishery Independent ---------------------
## setup -------------
# create spatial dataframe
hake_sdmTMB_df = hake_weight_age_df_updated %>% 
  dplyr::select(new_age, weight, catch_year, catch_month, cohort, length, distance_fished, sex_description, hb_longitude, hb_latitude) 
hake_sdmTMB_df$hb_longitude = -1 * hake_sdmTMB_df$hb_longitude
hake_sdmTMB_df = st_as_sf(hake_sdmTMB_df, coords = c("hb_longitude", "hb_latitude"))
st_crs(hake_sdmTMB_df) = 4269
hake_sdmTMB_df = st_transform(hake_sdmTMB_df, 32610)
coords =  st_coordinates(hake_sdmTMB_df)
coords = coords/1000 #convert to km
hake_sdmTMB_df = as.data.frame(cbind(hake_sdmTMB_df, coords, 
                                     hb_latitude = hake_weight_age_df_updated$hb_latitude,
                                     hb_longitude = -hake_weight_age_df_updated$hb_longitude))
# assign male and female to unassigned sexes. sampled randomly. Because these are earlier age classes, I assume a 50/50 chance of each - reasonable to assume that negligible sex-specific mortality has not yet occurred.

#hake_sdmTMB_df = hake_sdmTMB_df %>% 
#  mutate(sex_complete = sex_description) 
#hake_sdmTMB_df[hake_sdmTMB_df$sex_complete == "Unknown/Not determined" ,]$sex_complete <- #sample(c("Male", "Female"), length(hake_sdmTMB_df[hake_sdmTMB_df$sex_complete == "Unknown#/Not determined" ,]$sex_complete), prob = c(0.5, 0.5), replace = TRUE)


# Set data class as appropriate
#hake_sdmTMB_df$catch_month = as.factor(hake_sdmTMB_df$catch_month)
hake_sdmTMB_df$catch_year = as.integer(as.character(hake_sdmTMB_df$catch_year))
hake_sdmTMB_df$cohort = as.integer(hake_sdmTMB_df$cohort)
#hake_sdmTMB_df$sex_complete = as.factor(hake_sdmTMB_df$sex_complete)
#hake_sdmTMB_df$sex_description = as.factor(hake_sdmTMB_df$sex_description)

# create mesh
mesh <- make_mesh(hake_sdmTMB_df, xy_cols = c("X", "Y"), cutoff = 40)
plot(mesh)
# ggplot() +
#   inlabru::gg(mesh$mesh) +
#   geom_point(data = hake_sdmTMB_df, aes(x = X, y = Y, col = n)) +
#   coord_equal()

# visualize weight data
jpeg(filename = "plots/data_exploration/weight_data_hist.jpeg", units="in", width=6, height=3, res = 300)
par(mfrow = c(1,2))
hist(hake_sdmTMB_df$weight, main = "hist. of weight")
hist(log(hake_sdmTMB_df$weight), main = "hist of log weight")
dev.off()

# visualize sex differences
jpeg(filename = "plots/data_exploration/weight-at-age_by_sex.jpeg", units="in", width=6, height=3, res = 300)
ggplot(hake_sdmTMB_df, aes(x = new_age, y = weight, col = sex_description)) +
  geom_jitter(alpha = 0.5)
dev.off()



## Pairwise covariate model exploration-----------------------------------

# For this model exploration, I will explore 3 nested models
# 1. m_age <- weight ~ age
# 2. m_age_cohort <- weight ~ age + cohort
# 3. m_age_year <- weight ~ age + year

# 1.
m_age = sdmTMB( 
  data = hake_sdmTMB_df,
  formula = weight ~ 1 + s(new_age),
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off"
)

# 2.
# smoother on cohort
m_age_cohort = sdmTMB( 
  data = hake_sdmTMB_df,
  formula = weight ~ 1 + s(new_age) + s(cohort),
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off"
)

# cohort as random effect
hake_sdmTMB_df$cohort = as.factor(hake_sdmTMB_df$cohort)
m_age_cohort_re = sdmTMB( 
  data = hake_sdmTMB_df,
  formula = weight ~ s(new_age) + (1 | cohort),
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997,
                 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014,
                 2016)
)

# plot the cohort random effects
varying_intercepts = tidy(m_age_cohort_re, "ranef", conf.int = TRUE)[,1:2]
varying_intercepts$term <- c(1973:2016)
jpeg(filename = "plots/nested_models/cohort_random_effects.jpeg", units="in", width=5, height=3, res = 300)
ggplot(varying_intercepts, aes(x = term, y = estimate)) +
  geom_point() +
  theme_classic() +
  labs(title = "cohort random effects")
dev.off()

# plot the smoother on cohort
jpeg(filename = "plots/nested_models/cohort_smoother.jpeg", units="in", width=5, height=3, res = 300)
plot_smooth(m_age_cohort, select =  2)
dev.off()

# 3  

# year modeled as a random walk 
hake_sdmTMB_df$catch_year <- as.integer(as.character(hake_sdmTMB_df$catch_year))
m_age_year_rw = sdmTMB( 
  data = hake_sdmTMB_df,
  formula = weight ~ 0 + s(new_age),
  mesh = mesh, 
  time_varying = ~ 1,
  time = "catch_year",
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997,
                 1999,2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016
  )
)


temp_info = as.data.frame(cbind(year = c(1986, 1989, 1992, 1995, 1998, 2001, 2003, 2005, 2007, 2009, 2011, 2012, 2013, 2015, 2017), temp_anomaly = c("NA", "NA", "NA", "Avg", "Hot", "Cold", "Cold", "Avg", "Cold", "Cold", "Avg","Cold", "Cold", "Hot", "Hot")))
temp_info$year = as.numeric(temp_info$year)
extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997,
               1999,2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016
)
est <- as.list(m_age_year_rw$sd_report, "Est") # extract rw estimates
year_rw_est = as.data.frame(cbind(year = 1986:2017, est = est$b_rw_t)) # create dataframe
year_rw_est = year_rw_est %>% # determine whether data was interpolated
  mutate(interpolated = year %in% extra_time)
year_rw_est = left_join(year_rw_est, temp_info, by = join_by(year)) 
year_rw_est[is.na(year_rw_est)] <- "NA"

jpeg(filename = "plots/nested_models/year_rw.jpeg", units="in", width=5, height=3, res = 300)
ggplot(year_rw_est, aes(x = year, y = est)) +
  geom_line(col = "gray") +
  geom_point(aes(shape = interpolated, col = temp_anomaly)) +
  scale_color_manual(values = c("black", "blue", "red", "dark gray")) +
  scale_shape_manual(values = c(19,4), guide = "none") +
  labs(title = "Year as Random Walk", y = "coefficient") +
  theme_classic()
dev.off()



# year modeled as a random effect/intercept
hake_sdmTMB_df$catch_year <- as.factor(hake_sdmTMB_df$catch_year)
m_age_year_re = sdmTMB( 
  data = hake_sdmTMB_df,
  formula = weight ~ s(new_age) + (1 | catch_year),
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  control = sdmTMBcontrol(newton_loops = 1)
)

varying_intercepts = tidy(m_age_year_re, "ranef", conf.int = TRUE)[,1:2]
varying_intercepts$term <- c(1986, 1989, 1992, 1995, 1998, 2001, 2003, 2005, 2007, 2009, 2011, 2012, 2013, 2015, 2017)
varying_intercepts$temp_anomaly = c("NA", "NA", "NA", "Avg", "Hot", "Cold", "Cold", "Avg", "Cold", "Cold", "Avg","Cold", "Cold", "Hot", "Hot")

jpeg(filename = "plots/nested_models/year_re_tempColor.jpeg", units="in", width=5, height=3, res = 300)
ggplot(varying_intercepts, aes(x = term, y = estimate)) +
  geom_line(col = "light gray") +
  geom_point(aes(col = temp_anomaly)) +
  scale_color_manual(values = c("black", "blue", "red", "gray")) +
  theme_classic() +
  labs(title = "Year as Random Effect", x = "year", y = "coefficient")
dev.off()

AIC(m_age, m_age_cohort, m_age_cohort_re, m_age_year_re, m_age_year_rw) # model with age and cohort as random effect was the best fit


## Best fit cohortRE ---------------------

m_age_cohort_re_s = sdmTMB( 
  data = hake_sdmTMB_df,
  formula = weight ~ s(new_age) + (1 | cohort),
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on",
  spatiotemporal = "off",
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016)
)

m_age_cohort_re_st = sdmTMB( 
  data = hake_sdmTMB_df,
  formula = weight ~ s(new_age) + (1 | cohort),
  mesh = mesh, 
  family = lognormal(link = "log"),
  time = "catch_year",
  spatial = "off",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016)
)

m_age_cohort_re_s_st = sdmTMB( 
  data = hake_sdmTMB_df,
  formula = weight ~ s(new_age) + (1 | cohort),
  mesh = mesh, 
  family = lognormal(link = "log"),
  time = "catch_year",
  spatial = "on",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997,
                 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014,
                 2016)
)

AIC(m_age_cohort_re_s_st, m_age_cohort_re_st, m_age_cohort_re_s, m_age_cohort_re)

#save(m_age_cohort_re_s_st, file = "results/sdmTMB/cohortRE_s_st_model.RData")

## Cohort smoother -------------
hake_sdmTMB_df$cohort = as.integer(as.character(hake_sdmTMB_df$cohort))
m_age_cohort_sm_s = sdmTMB( 
  data = hake_sdmTMB_df,
  formula = weight ~ s(new_age) + s(cohort),
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on",
  spatiotemporal = "off",
  control = sdmTMBcontrol(newton_loops = 1),
  #extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016)
  )

m_age_cohort_sm_st = sdmTMB( 
  data = hake_sdmTMB_df,
  formula = weight ~ s(new_age) + s(cohort),
  mesh = mesh, 
  family = lognormal(link = "log"),
  time = "catch_year",
  spatial = "off",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016))

m_age_cohort_sm_s_st = sdmTMB( 
  data = hake_sdmTMB_df,
  formula = weight ~ s(new_age) + s(cohort),
  mesh = mesh, 
  family = lognormal(link = "log"),
  time = "catch_year",
  spatial = "on",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016))

#save(m_age_cohort_sm_s_st, file = "results/sdmTMB/cohortsmoother_s_st_model.RData")



## Month and Sex effects -------

hake_sdmTMB_df_complete_sex = hake_sdmTMB_df[hake_sdmTMB_df$sex_description %in% c("Female", "Male"),]
hake_sdmTMB_df_complete_sex$sex_description = as.factor(hake_sdmTMB_df_complete_sex$sex_description)
mesh <- make_mesh(hake_sdmTMB_df_complete_sex, xy_cols = c("X", "Y"), cutoff = 40)

m_age_month_cohort = sdmTMB( 
  data = hake_sdmTMB_df_complete_sex,
  formula = weight ~ s(new_age) + s(cohort) + catch_month,
  mesh = mesh, 
  family = lognormal(link = "log"),
  time = "catch_year",
  spatial = "on",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016))

#save(m_age_month_cohort, file = "results/sdmTMB/m_age_month_cohort.RData")

m_age_month_cohort_sex = sdmTMB( 
  data = hake_sdmTMB_df_complete_sex,
  formula = weight ~ s(new_age) + s(cohort) + catch_month + sex_description,
  mesh = mesh, 
  family = lognormal(link = "log"),
  time = "catch_year",
  spatial = "on",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016, 2018, 2020))

hake_sdmTMB_df$sex_complete = as.factor(hake_sdmTMB_df$sex_complete)
m_age_month_cohort_sex_insmooth = sdmTMB( 
  data = hake_sdmTMB_df_complete_sex,
  formula = weight ~ s(new_age, by = sex_description) + s(cohort) + catch_month + sex_description,
  mesh = mesh, 
  family = lognormal(link = "log"),
  time = "catch_year",
  spatial = "off",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016, 2018, 2020))
#save(m_age_month_cohort_sex_insmooth, file = "results/sdmTMB/m_age_month_cohort_sex_insmooth.RData")

#save(m_age_month_cohort_sex, file = "results/sdmTMB/m_age_month_cohort_sex.RData")

m_age_month_cohort_sex_sp_only = sdmTMB( 
  data = hake_sdmTMB_df_complete_sex,
  formula = weight ~ s(new_age) + s(cohort) + catch_month + sex_description,
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016, 2018, 2020))

m_age_month_cohort_sex_no_sp_st = sdmTMB( 
  data = hake_sdmTMB_df_complete_sex,
  formula = weight ~ s(new_age) + s(cohort) + catch_month + sex_description,
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "off",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016, 2018, 2020))

visreg::visreg(m_age_month_cohort_sex_no_sp_st, data = m_age_month_cohort_sex_no_sp_st$data, scale = "response")

m_age_month_cohort_sex_st_only = sdmTMB( 
  data = hake_sdmTMB_df_complete_sex,
  formula = weight ~ s(new_age) + s(cohort) + catch_month + sex_description,
  mesh = mesh, 
  time = "catch_year",
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016, 2018, 2020))
#save(m_age_month_cohort_sex_st_only, file = "results/sdmTMB/m_age_month_cohort_sex_st_only.RData")




# Fishery - Dependent ------------------------------------

load("results/can_df.RData"); load("results/us_df_full.RData")
# merge btoh canada and us info
can_df = can_df %>% 
  mutate(country = "canada")
us_df_full = us_df_full %>% 
  mutate(country = "usa")
fishery_df = rbind(can_df,us_df_full)
fishery_df$fcatch_year = as.factor(fishery_df$catch_year)
fishery_df$sex_description <- as.factor(fishery_df$sex_description)
fishery_df$cohort <- as.integer(fishery_df$cohort)
fishery_df$country <- as.factor(fishery_df$country)

## Data exploration ---------
fishery_df %>% 
  group_by(catch_month) %>% 
  summarise(avg = mean(weight), sd = sd(weight)) %>% 
  ggplot(aes(x = catch_month, y = avg)) +
  geom_point() +
  #geom_errorbar(aes(ymax = avg + sd, ymin = avg - sd, col = as.factor(country))) +
  theme_classic()

# higher catches in the summer when the fish are more northward (maybe model as a GAM?).

View(fishery_df %>% 
  filter(catch_month %in% c(12,1,2,3)) %>% 
  group_by(Source, country, catch_month, catch_year) %>% 
  summarise(n = n()))

jpeg("plots/data_exploration/fishery_data/samplesize_by_country.jpeg", units="in", width=12, height=7, res=600)
ggplot(fishery_df, aes(x = catch_month, group = country)) +
  geom_bar(aes(fill = country)) +
  theme_classic()
dev.off()

jpeg("plots/data_exploration/fishery_data/weight-at-age_by_country.jpeg", units="in", width=7, height=4, res=300)
fishery_df %>% 
  group_by(country, new_age) %>% 
  summarise(avg = mean(weight), n = n()) %>% 
  ggplot(aes(x = new_age, y = avg, group = country)) +
  geom_point(aes(col = country, size = n)) +
  theme_classic()
dev.off()

# larger weight-at-age for fishes in Canada.

jpeg("plots/data_exploration/fishery_data/avg_weight-at-age_timeseries.jpeg", units="in", width=7, height=4, res=300)
fishery_df %>%  # timeseries of average weight for each age 
  group_by(new_age, catch_year) %>% 
  summarise(avg = mean(weight), n = n()) %>% 
  ggplot(aes(x = catch_year, y = avg, group = new_age)) +
  geom_line(aes(col = new_age))
dev.off()

jpeg("plots/data_exploration/fishery_data/cohort_weight-at-age.jpeg", units="in", width=7, height=4, res=300)
fishery_df %>%  # cohort effects
  group_by(cohort, new_age) %>% 
  summarise(avg = mean(weight), n = n()) %>% 
  ggplot(aes(x = cohort, y = avg, group = new_age)) +
  geom_line(aes(col = new_age)) +
  theme_classic()
dev.off()

jpeg("plots/data_exploration/fishery_data/weightXage_by_sex.jpeg", units="in", width=7, height=4, res=300)
ggplot(fishery_df, aes(x = new_age, y = weight)) +
  geom_jitter(aes(col = sex_description), alpha = 0.2) +
  theme_classic()
dev.off()

# so much variability in weight-at-age, especially in females. Will definitely include sex as a covariate, and should probably try as a by factor in gam (need to include as a main effect) as well as a linear predictor like in survey model

## Model building ----------
# age as a GAM
# cohort as a GAM
# sex 1) as a by factor in age GAM and 2) as a linear predictor
# month as a GAM


# I started off with just age as a predictor and it did not converge. Added sex as a by factor in age gam
fm_age_sex = sdmTMB( 
    data = fishery_df,
    formula = weight ~ s(new_age, by = sex_description) + sex_description,
    family = lognormal(link = "log"),
    spatial = "off",
    spatiotemporal = "off",
    control = sdmTMBcontrol(newton_loops = 1)
  )

# also tried sex as a linear predictor
fm_age_sex_linear = sdmTMB( 
  data = fishery_df,
  formula = weight ~ s(new_age) + sex_description,
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  control = sdmTMBcontrol(newton_loops = 1)
)

# add cohort
fm_age_sex_cohort = sdmTMB( 
  data = fishery_df,
  formula = weight ~ s(new_age, by = sex_description) + s(cohort) + sex_description,
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  control = sdmTMBcontrol(newton_loops = 1)
)

# add month
fm_age_sex_cohort_month = sdmTMB( 
  data = fishery_df,
  formula = weight ~ s(new_age, by = sex_description) + s(cohort) + sex_description + s(catch_month),
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  control = sdmTMBcontrol(newton_loops = 1)
)

# add year as a random effect
fm_age_sex_cohort_month_yearre = sdmTMB( 
  data = fishery_df,
  formula = weight ~ s(new_age, by = sex_description) + s(cohort) + sex_description + s(catch_month) + (1 | fcatch_year),
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  control = sdmTMBcontrol(newton_loops = 1)
)

fm_age_sex_month_yearre = sdmTMB( 
  data = fishery_df,
  formula = weight ~ s(new_age, by = sex_description) + sex_description + s(catch_month) + (1 | fcatch_year),
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  control = sdmTMBcontrol(newton_loops = 1)
) # 20,000 AIC values less than the one with cohort

# I tried adding country on top but it did not work - the undefined sex factor is characteristic of the entire canadian fishery, and so I think that is already capturing the difference seen in Canada

# the model with cohort, year, age, sex, and month is the best model by FAR (at the very least, 10k AIC values)

#  save(fm_age_sex_cohort_month_yearre, file = "results/sdmTMB/fm_age_sex_cohort_month_yearre.RData")



# data exploration --------------
fishery_year = fishery_df %>% 
  group_by(catch_month, country) %>% 
  summarise(n = n()) %>% 
  mutate(source = "fishery")

survey_year = hake_sdmTMB_df %>% 
  group_by(catch_month) %>% 
  summarise(n = n()) %>% 
  mutate(source = "survey")

data_year = as.data.frame(rbind(fishery_year, survey_year))
jpeg(filename = "plots/data_exploration/data_availability_month.jpeg", units="in", width=8, height=3, res = 300) 
ggplot(data_year, aes(x = as.factor(catch_month), y = source)) +
  geom_point(aes(size = n, fill = source), shape = 21, col = "black") +
  theme_classic() +
  labs(x = "Month", y  = "", size = "sample size") +
  scale_fill_manual(values=c("#6e7cb9", "#d2848d")) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12))
dev.off()
