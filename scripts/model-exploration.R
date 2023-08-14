load("results/hake_weight_age_df.RData")
library(dplyr)
library(sf)
library(sp)
library(rgdal)
library(ggplot2)
library(sdmTMB)
library(visreg)

# setup -------------
# create spatial dataframe
hake_sdmTMB_df = hake_weight_age_df %>% 
  dplyr::select(new_age, weight, catch_year, catch_month, cohort, length, distance_fished, sex_description,
                hb_longitude, hb_latitude) 
hake_sdmTMB_df$hb_longitude = -1 * hake_sdmTMB_df$hb_longitude
hake_sdmTMB_df = st_as_sf(hake_sdmTMB_df, coords = c("hb_longitude", "hb_latitude"))
st_crs(hake_sdmTMB_df) = 4269
hake_sdmTMB_df = st_transform(hake_sdmTMB_df, 32610)
coords =  st_coordinates(hake_sdmTMB_df)
coords = coords/1000
hake_sdmTMB_df = as.data.frame(cbind(hake_sdmTMB_df, coords, 
                                     hb_latitude = hake_weight_age_df$hb_latitude,
                                     hb_longitude = hake_weight_age_df$hb_longitude))

# Set data class as appropriate
#hake_sdmTMB_df$catch_month = as.factor(hake_sdmTMB_df$catch_month)
hake_sdmTMB_df$catch_year = as.integer(as.character(hake_sdmTMB_df$catch_year))
hake_sdmTMB_df$cohort = as.factor(hake_sdmTMB_df$cohort)

# visualize weight data
jpeg(filename = "plots/data_exploration/weight_data_hist.jpeg", units="in", width=6, height=3, res = 300)
par(mfrow = c(1,2))
hist(hake_sdmTMB_df$weight, main = "hist. of weight")
hist(log(hake_sdmTMB_df$weight), main = "hist of log weight")
dev.off()

# create mesh
mesh <- make_mesh(hake_sdmTMB_df, xy_cols = c("X", "Y"), cutoff = 10)
plot(mesh)
# ggplot() +
#   inlabru::gg(mesh$mesh) +
#   geom_point(data = hake_sdmTMB_df, aes(x = X, y = Y, col = n)) +
#   coord_equal()

# Pairwise covariate model exploration-----------------------------------

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


# Best fit ---------------------

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

