load("results/hake_weight_age_df.RData")
load("results/sdmTMB/m4.alt.st.samplingyrs.RData")
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

# continuation of models
#hake_sdmTMB_df$catch_month = as.factor(hake_sdmTMB_df$catch_month)
hake_sdmTMB_df$catch_year = as.integer(as.character(hake_sdmTMB_df$catch_year))
hake_sdmTMB_df$cohort = as.integer(hake_sdmTMB_df$cohort)

# create mesh
mesh <- make_mesh(hake_sdmTMB_df, xy_cols = c("X", "Y"), cutoff = 10)
plot(mesh)
# ggplot() +
#   inlabru::gg(mesh$mesh) +
#   geom_point(data = hake_sdmTMB_df, aes(x = X, y = Y, col = n)) +
#   coord_equal()


# try diff models --------
hake_sdmTMB_df$cohort = as.integer(as.character(hake_sdmTMB_df$cohort))
hake_sdmTMB_df$catch_year = as.integer(as.character(hake_sdmTMB_df$catch_year))

m1 <- sdmTMB(
  data = hake_sdmTMB_df,
  formula = weight ~ s(new_age) + s(cohort),
  mesh = mesh, # can be omitted for a non-spatial model
  family = lognormal(link = "log"),
  spatial = "off"
)

m2 <- sdmTMB(
  data = hake_sdmTMB_df,
  formula = weight ~ s(new_age) + s(cohort),
  mesh = mesh, # can be omitted for a non-spatial model
  family = lognormal(link = "log"),
  spatial = "on",
  time = "catch_year",
  spatiotemporal = "AR1"
)
#save(m2, file = "results/sdmTMB/m2.RData")

m3 <- sdmTMB(
  data = hake_sdmTMB_df,
  formula = weight ~ s(new_age) + s(cohort) + catch_month,
  mesh = mesh, # can be omitted for a non-spatial model
  family = lognormal(link = "log"),
  spatial = "on",
  time = "catch_year",
  spatiotemporal = "AR1"
)
#save(m3, file = "results/sdmTMB/m3.RData")

m4 <- sdmTMB(
  data = hake_sdmTMB_df,
  formula = weight ~ 0 + s(new_age) + s(cohort) + catch_month,
  time_varying = ~ 1,
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on",
  time = "catch_year",
  spatiotemporal = "AR1"
)

m5 <- sdmTMB(
  data = hake_sdmTMB_df,
  formula = weight ~ 0 + s(new_age) + s(cohort) + catch_month + s(catch_year),
  time_varying = ~ 1,
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on",
  time = "catch_year",
  spatiotemporal = "AR1"
)

hake_sdmTMB_df$catch_year = as.factor(hake_sdmTMB_df$catch_year)
m6 <- sdmTMB(
  data = hake_sdmTMB_df,
  formula = weight ~ 0 + s(new_age, by = catch_year)  + s(cohort) + catch_month,
  time_varying = ~ 1,
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on",
  time = "catch_year",
  spatiotemporal = "AR1"
)

m7 <- sdmTMB(
  data = hake_sdmTMB_df,
  formula = weight ~ 0 + s(new_age, by = catch_year) + catch_month,
  time_varying = ~ 1,
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on",
  time = "catch_year",
  spatiotemporal = "AR1"
)

m8 <- sdmTMB(
  data = hake_sdmTMB_df,
  formula = weight ~ new_age,
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on"
)



# model diagnostics
tidy(m8, "ran_pars", conf.int = TRUE)
hake_sdmTMB_df$resids_m8_simple <- residuals(m8) # randomized quantile residuals
predict_m8 = predict(m8)
#set.seed(1)
#rq_res <- residuals(m8) # randomized quantile residuals
#qqnorm(rq_res);qqline(rq_res)

#jpeg(file="plots/sdmTMB/residuals_m8_simple.jpeg")
ggplot(hake_sdmTMB_df, aes(X, Y, col = resids_m8_simple)) +
  scale_colour_gradient2() +
  geom_point() +
  facet_wrap(~catch_year) +
  coord_fixed() +
  labs(title= "weight ~ new_age as spatial model") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
#dev.off()


est <- as.list(m5$sd_report, "Est")
plot(est$b_rw_t)





# random plotting ---------

hake_sdmTMB_df %>% 
  group_by(new_age, catch_year) %>% 
  summarise(avg_weight = mean(weight)) %>% 
  ggplot(aes(x = catch_year, y = avg_weight, group = as.factor(new_age), col = as.factor(new_age))) +
      geom_line() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))

hake_sdm_north = subset(hake_sdmTMB_df, Y > 4957.963)
hake_sdm_south = subset(hake_sdmTMB_df, Y < 4957.963)
p1 = hake_sdm_north %>% 
  group_by(new_age, catch_year) %>% 
  summarise(avg_weight = mean(weight)) %>% 
  ggplot(aes(x = catch_year, y = avg_weight, group = as.factor(new_age), col = as.factor(new_age))) +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "North") +
  lims(y = c(0,4))
p2 = hake_sdm_south %>% 
  group_by(new_age, catch_year) %>% 
  summarise(avg_weight = mean(weight)) %>% 
  ggplot(aes(x = catch_year, y = avg_weight, group = as.factor(new_age), col = as.factor(new_age))) +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "South") +
  lims(y = c(0,4))

ggarrange(p1,p2,labels = c("A", "B"), common.legend = TRUE, legend = "right")

#jpeg(file="plots/weight-at-age/weight_at_age_by_year_north.jpeg")
hake_sdmTMB_df %>% 
  filter( Y < 4957.963) %>% 
  group_by(new_age, catch_year) %>% 
  summarise(avg_weight = mean(weight), n =n()) %>% 
  ggplot(aes(x = new_age, y = avg_weight, group = as.factor(catch_year), col = as.factor(catch_year))) +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
#dev.off()

#jpeg(file="plots/weight-at-age/observations_of_age_by_month.jpeg")
hake_sdmTMB_df %>% 
  group_by(catch_month, new_age) %>% 
  summarise(avg_weight = mean(weight), n = n()) %>% 
  ggplot(aes(x = catch_month, y = n, group = as.factor(new_age), fill = as.factor(new_age))) +
  geom_col() +
  theme_classic() 
#dev.off()




# continuation of models
hake_sdmTMB_df$catch_month = as.factor(hake_sdmTMB_df$catch_month)
hake_sdmTMB_df$catch_year = as.integer(as.character(hake_sdmTMB_df$catch_year))
hake_sdmTMB_df$cohort = as.integer(hake_sdmTMB_df$cohort)

m4.spatiotemporal = sdmTMB(
  data = hake_sdmTMB_df,
  formula = weight ~ 0 + s(new_age) + s(cohort) + (1 | catch_month) + (1 | catch_year),
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on",
  time = "catch_year",
  spatiotemporal = "AR1"
)
#save(m4.spatiotemporal, file = "results/sdmTMB/m4.spatiotemporal.RData")

m4.alt.spatiotemporal = sdmTMB(
  data = hake_sdmTMB_df,
  formula = 0 + weight ~ s(new_age) + s(cohort) + (1 | catch_month),
  time_varying = ~1,
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002,
                 2004, 2006, 2008, 2010, 2014, 2016),
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on",
  time = "catch_year",
  spatiotemporal = "AR1"
) # Can't predict missing or future years if catch_year is a factor. Instead 
# switched to a time-varying intercept. 

# This it the one!!! ----------------
m4.alt.st.samplingyrs = sdmTMB(
  data = hake_sdmTMB_df,
  formula = weight ~ 0 + s(new_age) + s(cohort) + (1 | catch_month),
  mesh = mesh, 
  time_varying = ~ 1,
  family = lognormal(link = "log"),
  spatial = "on",
  time = "catch_year",
  spatiotemporal = "AR1"
)
#save(m4.alt.st.samplingyrs, file = "results/sdmTMB/m4.alt.st.samplingyrs.RData")

m4.spatial.only = sdmTMB(
  data = hake_sdmTMB_df,
  formula = weight ~ 0 + s(new_age) + s(cohort) + (1 | catch_month) + (1 | catch_year),
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on"
)

m4.no.spatiotemp = sdmTMB(
  data = hake_sdmTMB_df,
  formula = weight ~ 0 +  s(new_age) + s(cohort) + (1 | catch_month) + (1 | catch_year),
  family = lognormal(link = "log"),
  spatial = "off"
)

coefficients_m4.alt.st.samplingyrs = tidy(m4.alt.st.samplingyrs, "ran_pars", conf.int = TRUE)
hake_sdmTMB_df$resids_m4spatiotemporal <- residuals(m4.spatiotemporal) # randomized quantile residuals
predict_m4.st.alt = predict(m4.alt.st.samplingyrs)
#save(predict_m4.spatiotemporal, file = "results/sdmTMB/predict_m4.spatiotemporal.RData")

coefficients_m4.alt.st.samplingyrs %>% 
  filter(!term %in% "range") %>% 
ggplot(aes( x = estimate, y = term)) +
  geom_point() +
  theme_classic() +
  lims( x = c(-0.1, 0.4)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high), width=.1) +
  labs(y = "coefficient", subtitle = "point estimates with confidence intervals")

# check residuals
rq_res <- residuals(m4.spatiotemporal)
rq_res <- rq_res[is.finite(rq_res)] # some Inf
qqnorm(rq_res);qqline(rq_res)

mcmc_res <- residuals(m4.spatiotemporal, type = "mle-mcmc", mcmc_iter = 201, mcmc_warmup = 200) 
qqnorm(mcmc_res);qqline(mcmc_res)

jpeg(filename = "plots/sdmTMB/m4.st.alt.spatialeffects.jpeg")
ggplot(predict_m4.st.alt, aes(X, Y, col = omega_s)) +
  scale_colour_gradient2() +
  geom_point() +
  #facet_wrap(~catch_year) +
  coord_fixed() +
  labs(title= "m4 spatiotemporal") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# making predictions - projecting to the whole region --------------------------------------

# makes all combinations of x and y:
load("results/grid_pred.Rdata")
load("results/sdmTMB/m4.2.RData")

year_vector = sort(as.numeric(unique(hake_weight_age_df$catch_year)))
month_vector = sort(as.numeric(unique(hake_weight_age_df$catch_month)))
age_vector = sort(as.numeric(unique(hake_weight_age_df$new_age)))
cohort_vector = sort(as.numeric(unique(hake_weight_age_df$cohort)))
grid_pred$catch_year = year_vector[1]
grid_pred_sdm = grid_pred
for(i in 2:length(year_vector)) { # add years
  grid_pred$catch_year = year_vector[i]
  grid_pred_sdm = rbind(grid_pred_sdm, grid_pred)
}
grid_pred_sdm$catch_month = month_vector[1] # add month
grid_pred_sdm_full = grid_pred_sdm
for(i in 2:length(month_vector)) { 
  grid_pred$catch_month = month_vector[i]
  grid_pred_sdm_full = rbind(grid_pred_sdm_full, grid_pred_sdm)
}
grid_pred_sdm_full$new_age = age_vector[1] # add age
grid_pred_sdm_full2 = grid_pred_sdm_full
for(i in 2:length(age_vector)) { 
  grid_pred_sdm_full$new_age = age_vector[i]
  grid_pred_sdm_full2 = rbind(grid_pred_sdm_full2, grid_pred_sdm_full)
}
grid_pred_sdm_full2$cohort = cohort_vector[1] # add age
grid_pred_sdm_full3 = grid_pred_sdm_full2
for(i in 2:length(cohort_vector)) { 
  grid_pred_sdm_full2$cohort = cohort_vector[i]
  grid_pred_sdm_full3 = rbind(grid_pred_sdm_full3, grid_pred_sdm_full2)
}


grid_pred_sdm_full3$catch_month = as.integer(as.character(grid_pred_sdm_full3$catch_month))
grid_pred_sdm_full3$new_age = as.integer(as.character(grid_pred_sdm_full3$new_age))
grid_pred_sdm_full3$X = grid_pred_sdm_full3$X/1000
grid_pred_sdm_full3$Y = grid_pred_sdm_full3$Y/1000
grid_pred_sdm_full3 = grid_pred_sdm_full3[,c(1:2,4:7)]
colnames(grid_pred_sdm_full3) = c("X", "Y", "catch_year", "catch_month", "new_age", "cohort")
plot(grid_pred_sdm_full3$X, grid_pred_sdm_full3$Y)

predicted_vals = predict(m4.5, newdata = grid_pred_sdm_full)



# Varying intercepts
temp_anomaly = c("NA", "NA", "NA", "Avg", "Hot", "Cold", "Cold", "Avg", 
                 "Cold", "Cold", "Avg","Cold", "Cold", "Hot", "Hot")
est <- as.list(m4.2$sd_report, "Est")
est_se <- as.list(m4.2$sd_report, "Std. Error")
varying_intercept = as.data.frame(cbind(intercept = est$b_rw_t, se = est_se$b_rw_t, year_vector, temp = temp_anomaly))
varying_intercept$intercept = as.numeric(varying_intercept$intercept)
varying_intercept$se = as.numeric(varying_intercept$se)

jpeg(filename = "plots/sdmTMB/time-varying-intercepts.m4.2.jpeg", units="in", width=5, height=5, res = 300)
ggplot(varying_intercept, aes(x = year_vector, y = intercept, col = temp)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values=c("grey", "blue", "red", "black")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Time-Varying Intercept - m4.2")
dev.off()

  

# Highlight everything and run it - otherwise won't run.
jpeg(filename = "plots/sdmTMB/month-varying-intercepts.jpeg", units="in", width=4, height=3, res = 300)
ggplot(month_vary_intercept, aes(x = term, y = estimate)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Month Time-Varying Intercept") + 
  scale_x_discrete(labels=c("catch_month_6" = "June", 
                            "catch_month_7" = "July",
                            "catch_month_8" = "August",
                            "catch_month_9" = "September")) +
  theme_classic()
dev.off()


# Visualising GAMs
jpeg(filename = "plots/sdmTMB/cohort-GAM-scaled.jpeg", units="in", width=4, height=3, res = 300)
visreg::visreg(m4.5, xvar = "cohort", 
               xlim = c(1980, 2015), ylim = c(-0.1, 2), scale = "response")
title(main = "cohort smoothed function")
dev.off()

jpeg(filename = "plots/sdmTMB/cohort-smoothed-function.jpeg", units="in", width=4, height=3, res = 300)
visreg::visreg(m4.5, xvar = "cohort", 
               xlim = c(1980, 2015), ylim = c(-5, 5))
title(main = "cohort smoothed function")
dev.off()

jpeg(filename = "plots/sdmTMB/weight-at-age-GAM.jpeg", units="in", width=4, height=3, res = 300)
visreg::visreg(m4.5, xvar = "new_age", 
               xlim = c(0, 15), ylim = c(-0.1, 2), scale = "response")
title(main = "weight-at-age smoothed function - scaled")
dev.off()

jpeg(filename = "plots/sdmTMB/weight-at-age-GAM-notscaled.jpeg", units="in", width=4, height=3, res = 300)
visreg::visreg(m4.5, xvar = "new_age", xlim = c(0, 15), ylim = c(-7, 7))
title(main = "weight-at-age smoothed function")
dev.off()


# okay let's revisit the model and make catch_month a linear predictor.

m4.2 = sdmTMB( # spatial and spatiotemporal
  data = hake_sdmTMB_df,
  formula = weight ~ 0 + s(new_age) + s(cohort) + catch_month,
  mesh = mesh, 
  time_varying = ~ 1,
  family = lognormal(link = "log"),
  spatial = "on",
  time = "catch_year",
  spatiotemporal = "AR1"
)
#save(m4.2, file = "results/sdmTMB/m4.2.RData")
pred_avg_sim <- predict(m4.2, return_tmb_object = TRUE)
index_avg_sim <- get_index(pred_avg_sim, bias_correct = TRUE)
cond_index <- ggplot(index_avg_sim, aes(x = catch_year, y = est)) +
  labs(y = "Predicted weight", x = "Year") +
  geom_point(size = 2.5, shape = 21, fill = "gray1", color = "white") +
  geom_errorbar(data = index_avg_sim, inherit.aes = FALSE, 
                aes(x = catch_year, ymax = upr, ymin = lwr),
                width = 0, alpha = 0.8, color = "gray1", size = 0.7) +
  theme(axis.title = element_text(size = 9),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.title = element_text(size = 9),
        legend.spacing.x = unit(0, 'cm'),
        plot.margin = unit(c(0.4, 0.4, 0.4, 0), "cm")) +
  guides(linetype = "none",
         color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1,
                              keywidth = 0.1, keyheight = 0.1, default.unit = "inch",
                              override.aes = list(linetype = c(1,1,1,1,1,0),
                                                  shape = c(NA,NA,NA,NA,NA,16),
                                                  alpha = rep(0.8, 6),
                                                  color = pal)))

m4.3 = sdmTMB( # just spatiotemporal no spatial
  data = hake_sdmTMB_df,
  formula = weight ~ 0 + s(new_age) + s(cohort) + catch_month,
  mesh = mesh, 
  time_varying = ~ 1,
  family = lognormal(link = "log"),
  spatial = "off",
  time = "catch_year",
  spatiotemporal = "AR1"
) # this has very similar AIC value to m4.2, but the sanity check revealed some 
# convergence issues

m4.4 = sdmTMB( # spatial and spatiotemporal
  data = hake_sdmTMB_df,
  formula = weight ~ 0 + s(new_age) + s(cohort) + catch_month + as.factor(catch_year),
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on",
  time = "catch_year",
  spatiotemporal = "AR1"
) # did not converge


hake_sdmTMB_df$catch_year = as.factor(hake_sdmTMB_df$catch_year)
m4.5 = sdmTMB( # spatial and spatiotemporal
  data = hake_sdmTMB_df,
  formula = weight ~ 1 + s(new_age) + s(cohort) + catch_month + (1 | catch_year),
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "on",
  time = "catch_year",
  spatiotemporal = "AR1",
  control = sdmTMBcontrol(newton_loops = 1)
)

temp_anomaly = c("NA", "NA", "NA", "Avg", "Hot", "Cold", "Cold", "Avg", 
                 "Cold", "Cold", "Avg","Cold", "Cold", "Hot", "Hot")
jpeg(filename = "plots/sdmTMB/m4.5-year_RE.jpeg", units="in", width=4, height=3, res = 300)
tidy(m4.5, "ranef", conf.int = TRUE) %>% 
  mutate(year = sort(unique(hake_sdmTMB_df$catch_year)), temp = temp_anomaly) %>% 
  ggplot(aes(x = year, y = estimate, col = temp)) +
  geom_point() +  
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2) +
  scale_color_manual(values=c("grey", "blue", "red", "black")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "catch_year as random effect")
dev.off()

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
  spatiotemporal = "off"
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

# 3  year modeled as a random effect/intercept
hake_sdmTMB_df$catch_year <- as.factor(hake_sdmTMB_df$catch_year)
m_age_year = sdmTMB( 
  data = hake_sdmTMB_df,
  formula = weight ~ s(new_age) + (1 | catch_year),
  mesh = mesh, 
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  control = sdmTMBcontrol(newton_loops = 1)
)

varying_intercepts = tidy(m_age_year, "ranef", conf.int = TRUE)[,1:2]
varying_intercepts$term <- c(1986, 1989, 1992, 1995, 1998, 2001, 2003, 2005, 2007, 2009, 2011, 2012, 2013, 2015, 2017)
varying_intercepts$temp_anomaly = c("NA", "NA", "NA", "Avg", "Hot", "Cold", "Cold", "Avg", "Cold", "Cold", "Avg","Cold", "Cold", "Hot", "Hot")

ggplot(varying_intercepts, aes(x = term, y = estimate, col = temp_anomaly)) +
  geom_point() +
  scale_color_manual(values = c("gray", "blue", "red", "black")) +
  theme_classic()