load("results/hake_weight_age_df.RData")
library(dplyr)
library(sf)
library(sp)
library(rgdal)
library(ggplot2)
library(sdmTMB)


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

# create mesh
mesh <- make_mesh(hake_sdmTMB_df, xy_cols = c("X", "Y"), cutoff = 10)
plot(mesh)
# ggplot() +
#   inlabru::gg(mesh$mesh) +
#   geom_point(data = hake_sdmTMB_df, aes(x = X, y = Y, col = n)) +
#   coord_equal()


# try diff models
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

coefficients_m4_spatiotemporal = tidy(m4.spatiotemporal, "ran_pars", conf.int = TRUE)
hake_sdmTMB_df$resids_m4spatiotemporal <- residuals(m4.spatiotemporal) # randomized quantile residuals
predict_m4.st.alt = predict(m4.alt.st.samplingyrs)
#save(predict_m4.spatiotemporal, file = "results/sdmTMB/predict_m4.spatiotemporal.RData")

coefficients_m4_spatiotemporal %>% 
  filter(!term %in% "range") %>% 
ggplot(aes( x = estimate, y = term)) +
  geom_point() +
  theme_classic() +
  lims( x = c(-0.25, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

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
load("results/sdmTMB/m4.spatiotemporal.RData")

year_vector = sort(as.numeric(unique(hake_weight_age_df$catch_year)))
grid_pred$catch_year = year_vector[1]
grid_pred_sdm = grid_pred
for(i in 2:length(year_vector)) {
  grid_pred$catch_year = year_vector[i]
  grid_pred_sdm = rbind(grid_pred_sdm, grid_pred)
}
grid_pred_sdm = grid_pred_sdm[,c(1,2,4)]
colnames(grid_pred_sdm) = c("X", "Y", "catch_year")
grid_pred_sdm$catch_year = as.integer(as.character(grid_pred_sdm$catch_year))
grid_pred_sdm$X = grid_pred_sdm$X/1000
grid_pred_sdm$Y = grid_pred_sdm$Y/1000

predicted_vals = predict(m4.alt.st.samplingyrs, newdata = grid_pred_sdm)



# Varying intercepts
temp_anomaly = c("NA", "NA", "NA", "Avg", "Hot", "Cold", "Cold", "Avg", 
                 "Cold", "Cold", "Avg","Cold", "Cold", "Hot", "Hot")
est <- as.list(m4.alt.st.samplingyrs$sd_report, "Est")
est_se <- as.list(m4.alt.st.samplingyrs$sd_report, "Std. Error")
varying_intercept = as.data.frame(cbind(intercept = est$b_rw_t, se = est_se$b_rw_t, year_vector, temp = temp_anomaly))
varying_intercept$intercept = as.numeric(varying_intercept$intercept)
varying_intercept$se = as.numeric(varying_intercept$se)

jpeg(filename = "plots/sdmTMB/time-varying-intercepts.jpeg")
ggplot(varying_intercept, aes(x = year_vector, y = intercept, col = temp)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values=c("grey", "blue", "red", "black")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Time-Varying Intercept")
dev.off()


