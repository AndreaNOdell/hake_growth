load("results/hake_weight_age_df.RData")
library(dplyr)
library(sf)
library(sp)
library(rgdal)
library(ggplot2)
library(sdmTMB)


# create spatial dataframe
hake_sdmTMB_df = hake_weight_age_df %>% 
  select(new_age, weight, catch_year, catch_month, cohort, length, distance_fished, sex_description,
         hb_longitude, hb_latitude) 
hake_sdmTMB_df$hb_longitude = -1 * hake_sdmTMB_df$hb_longitude
hake_sdmTMB_df = st_as_sf(hake_sdmTMB_df, coords = c("hb_longitude", "hb_latitude"))
st_crs(hake_sdmTMB_df) = 4269
hake_sdmTMB_df = st_transform(hake_sdmTMB_df, 32610)
coords =  st_coordinates(hake_sdmTMB_df)
coords = coords/1000
hake_sdmTMB_df = as.data.frame(cbind(hake_sdmTMB_df, coords))

# create mesh
mesh <- make_mesh(hake_sdmTMB_df, xy_cols = c("X", "Y"), cutoff = 1)
plot(mesh)
# ggplot() +
#   inlabru::gg(mesh$mesh) +
#   geom_point(data = hake_sdmTMB_df, aes(x = X, y = Y, col = n)) +
#   coord_equal()


# try diff models
hake_sdmTMB_df$cohort = as.numeric(hake_sdmTMB_df$cohort)

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

# model diagnostics
tidy(m3, "ran_pars", conf.int = TRUE)
hake_sdmTMB_df$resids_m3 <- residuals(m3)# randomized quantile residuals
hake_sdmTMB_df$predict_m3 = predict(m3)
#qqnorm(hake_sdmTMB_df$resids)
#qqline(hake_sdmTMB_df$resids)

#jpeg(file="plots/sdmTMB/spatiotemporal_RE_m3.jpeg")
ggplot(hake_sdmTMB_df, aes(X, Y, col = predict_m3$epsilon_st)) +
  scale_colour_gradient2() +
  geom_point() +
  facet_wrap(~catch_year) +
  coord_fixed() +
  labs(title= "spatiotemporal random effects")
#dev.off()







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


hake_sdmTMB_df %>% 
  group_by(cohort, catch_year) %>% 
  summarise(avg_weight = mean(weight)) %>% 
  ggplot(aes(x = catch_year, y = avg_weight, group = as.factor(cohort), col = as.factor(cohort))) +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


