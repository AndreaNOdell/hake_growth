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

# Survey -----------
# Load data 
load("results/sdmTMB/m_age_month_cohort_sex_insmooth.RData")
load("results/hake_weight_age_df_updated.RData")
model = m_age_month_cohort_sex_insmooth
dat = model$data

# get predictions for current observations
pred <- predict(model)

# Gather unique combinations of covariates
dat$year_age_cohort_month_sex <- as.factor(paste(dat$catch_year, dat$new_age, dat$cohort, dat$catch_month, dat$sex_description)) #, dat$sex_complete, dat$X, dat$Y

# Create spatial grid
grid_pred_full <- pred_grid(x = dat$X*1000, y = dat$Y*1000,
                            res =5e3, plot = FALSE)
names(grid_pred_full) <- c("easting_m", "northing_m")
grid_pred_full$cell_id <- 1:nrow(grid_pred_full)
#plot(grid_pred_full$easting_m, grid_pred_full$northing_m)


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
#plot(grid_pred$easting_m, grid_pred$northing_m)

grid_pred = grid_pred %>% 
  mutate(coords = paste(easting_m, northing_m)) %>% 
  select(!cell_id)

# Expand covariate combinations
df <- expand.grid(year_age_cohort_month_sex = unique(dat$year_age_cohort_month_sex))
df$catch_year <- as.numeric(unlist(lapply(strsplit(as.character(df$year_age_cohort_month_sex), " "), getElement, 1)))
df$new_age <- as.numeric(unlist(lapply(strsplit(as.character(df$year_age_cohort_month_sex), " "), getElement, 2)))
df$cohort <- unlist(lapply(strsplit(as.character(df$year_age_cohort_month_sex), " "), getElement, 3))
df$catch_month <- unlist(lapply(strsplit(as.character(df$year_age_cohort_month_sex), " "), getElement, 4))
df$sex_description <- unlist(lapply(strsplit(as.character(df$year_age_cohort_month_sex), " "), getElement, 5))
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
pred_grid$sex_description <- as.factor(pred_grid$sex_description)
pred_grid$catch_year <- as.numeric(pred_grid$catch_year)
#pred_grid$X <- as.numeric(pred_grid$X)
#pred_grid$Y <- as.numeric(pred_grid$Y)

pred_out <- predict(model, newdata = pred_grid, return_tmb_object = FALSE)#, se_fit = TRUE)
years_sampled = as.numeric(unique(hake_weight_age_df_updated$catch_year))
pred_out = pred_out[pred_out$catch_year %in% years_sampled,]
coords = cbind(pred_out$X, pred_out$Y)
sputm <- SpatialPoints(coords, proj4string=CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
spgeo <- as.data.frame(spTransform(sputm, CRS("+proj=longlat +datum=WGS84")))
coords = cbind(as.data.frame(sputm), spgeo)
names(coords) = c("X_test", "Y_test", "hb_longitude", "hb_latitude")
pred_out = cbind(pred_out, coords[,3:4]) 

# Fishery -------------
load("results/sdmTMB/fm_age_sex_cohort_month_yearre.RData")

model = fm_age_sex_cohort_month_yearre
dat = model$data

pred_out_fishery = predict(fm_age_sex_cohort_month_yearre)
save(pred_out_fishery, file = "results/sdmTMB/pred_out_fishery.RData")


