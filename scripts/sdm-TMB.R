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
  spatiotemporal = "IID"
)

# model diagnostics
tidy(m2, "ran_pars", conf.int = TRUE)
hake_sdmTMB_df$resids <- residuals(m1) # randomized quantile residuals
qqnorm(hake_sdmTMB_df$resids)
qqline(hake_sdmTMB_df$resids)

jpeg(file="plots/sdmTMB/spatiotemp_cohort_resids_by_year.jpeg")
ggplot(hake_sdmTMB_df, aes(X, Y, col = resids)) +
  scale_colour_gradient2() +
  geom_point() +
  facet_wrap(~catch_year) +
  coord_fixed()
dev.off()
