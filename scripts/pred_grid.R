library(dplyr)
library(sf)
library(sp)
library(rgdal)
library(ggplot2)
library(sdmTMB)
library(visreg)
library(ggpubr)
library(ggbreak)
source("malick_code/functions.R")

# Load data
load("results/sdmTMB/m_age_month_cohort.RData")
load("results/hake_weight_age_df.RData")
model = m_age_month_cohort
dat = model$data

# get predictions for current observations
pred <- predict(model)

# Gather unique combinations of covariates
dat$year_age_cohort_month_sex <- as.factor(paste(dat$catch_year, dat$new_age, dat$cohort, dat$catch_month)) #, dat$sex_complete, dat$X, dat$Y

# Create spatial grid
grid_pred_full <- pred_grid(x = dat$X*1000, y = dat$Y*1000,
                            res =5e3, plot = FALSE)
names(grid_pred_full) <- c("easting_m", "northing_m")
grid_pred_full$cell_id <- 1:nrow(grid_pred_full)
plot(grid_pred_full$easting_m, grid_pred_full$northing_m)


excl <- exclude_grid(grid_pred_full$easting_m,
                     grid_pred_full$northing_m,
                     dat$X*1000,
                     dat$Y*1000,
                     method = "euclidean",
                     dist = 10e3,
                     plot = FALSE)

grid_pred <- grid_pred_full[!excl, ]
row.names(grid_pred) <- NULL

## Plot grid
plot(grid_pred$easting_m, grid_pred$northing_m)

grid_pred = grid_pred %>% 
  mutate(coords = paste(easting_m, northing_m)) %>% 
  select(!cell_id)

# Expand covariate combinations
df <- expand.grid(year_age_cohort_month_sex = unique(dat$year_age_cohort_month_sex))
df$catch_year <- as.numeric(unlist(lapply(strsplit(as.character(df$year_age_cohort_month_sex), " "), getElement, 1)))
df$new_age <- as.numeric(unlist(lapply(strsplit(as.character(df$year_age_cohort_month_sex), " "), getElement, 2)))
df$cohort <- unlist(lapply(strsplit(as.character(df$year_age_cohort_month_sex), " "), getElement, 3))
df$catch_month <- unlist(lapply(strsplit(as.character(df$year_age_cohort_month_sex), " "), getElement, 4))
#df$X <- unlist(lapply(strsplit(as.character(df$year_age_cohort_month_sex), " "), getElement, 5))
#df$Y <- unlist(lapply(strsplit(as.character(df$year_age_cohort_month_sex), " "), getElement, 6))
df$cell = 1:nrow(df)

# now expand the spatial grid so there's all the spatial info for each cell of the df

spat_df = expand.grid(cell = df$cell, coords = grid_pred$coords) 
spat_df$X =  as.numeric(unlist(lapply(strsplit(as.character(spat_df$coords), " "), getElement, 1)))/1000
spat_df$Y =  as.numeric(unlist(lapply(strsplit(as.character(spat_df$coords), " "), getElement, 2)))/1000

# Then join covariate grid with spatial grid
pred_grid = full_join(spat_df, df)

pred_grid$new_age <- as.integer(pred_grid$new_age)
pred_grid$cohort <- as.integer(pred_grid$cohort)
pred_grid$catch_month <- as.numeric(as.character(pred_grid$catch_month))
#pred_grid$sex_complete <- as.factor(pred_grid$sex_complete)
pred_grid$catch_year <- as.numeric(pred_grid$catch_year)
#pred_grid$X <- as.numeric(pred_grid$X)
#pred_grid$Y <- as.numeric(pred_grid$Y)

pred_out <- predict(model, newdata = pred_grid, return_tmb_object = FALSE)#, se_fit = TRUE)
coords = cbind(pred_out$X, pred_out$Y)
sputm <- SpatialPoints(coords, proj4string=CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
spgeo <- as.data.frame(spTransform(sputm, CRS("+proj=longlat +datum=WGS84")))
coords = cbind(as.data.frame(sputm), spgeo)
names(coords) = c("X_test", "Y_test", "hb_longitude", "hb_latitude")
pred_out = cbind(pred_out, coords[,3:4]) 

# Index  ----------------
# Whole spatial domain
years_sampled = sort(unique(hake_weight_age_df$catch_year))
pred_sims <- predict(model, newdata = pred_grid, return_tmb_object = TRUE)
index_sims = get_index(pred_sims)
index_sims = index_sims %>% 
  filter(catch_year %in% years_sampled) %>% 
  mutate(region = "All")

n_peryear = pred_grid %>% 
  filter(catch_year %in% years_sampled) %>% 
  group_by(catch_year) %>% 
  summarise(n = n()) %>% 
  mutate(region = "All")


#  Canada
pred_grid_north =  pred_grid %>% 
  filter(Y > 5437.8)  # 547 is the UTM coordinate associated with 49 latitude
pred_sims_north = predict(model, newdata = pred_grid_north, return_tmb_object = TRUE)
index_sims_north = get_index(pred_sims_north)
index_sims_north = index_sims_north %>% 
  filter(catch_year %in% years_sampled) %>% 
  mutate(region = "Canada")

n_peryear_north = pred_grid_north %>% 
  filter(catch_year %in% years_sampled) %>% 
  group_by(catch_year) %>% 
  summarise(n = n()) %>% 
  mutate(region = "Canada")

# USA
pred_grid_south = pred_grid %>% 
  filter(Y < 5437.8) 
pred_sims_south = predict(model, newdata = pred_grid_south, return_tmb_object = TRUE)
index_sims_south = get_index(pred_sims_south) 
index_sims_south = index_sims_south %>% 
  filter(catch_year %in% years_sampled) %>% 
  mutate(region = "USA")

n_peryear_south = pred_grid_south %>% 
  filter(catch_year %in% years_sampled) %>% 
  group_by(catch_year) %>% 
  summarise(n = n()) %>% 
  mutate(region = "USA")

index_sims_full = rbind(index_sims, index_sims_north, index_sims_south)
n_peryear_full = rbind(n_peryear, n_peryear_north, n_peryear_south)
index_final = left_join(index_sims_full, n_peryear_full) %>% 
  mutate(est_weight = est/n, lwr_weight = lwr/n, upr_weight = upr/n)

#save(index_final, file = "results/sdmTMB/index_final-m_age_month_cohort.RData")


# plot the weight index by region
jpeg(filename = "plots/output_exploration/weight_index.jpeg", units="in", width=6, height=4, res = 300)
ggplot(index_final, aes(x = catch_year, y = est_weight)) +
  geom_ribbon(aes(ymin = lwr_weight, ymax = upr_weight, fill = region), alpha = 0.4) +
  geom_line(aes(color = region)) +
  theme_classic()
dev.off()








