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
dat$year_age_cohort_month_sex <- as.factor(paste(dat$catch_year, dat$new_age, dat$cohort, dat$catch_month)) #, dat$sex_complete

# Create spatial grid
grid_pred_full <- pred_grid(x = dat$X*1000, y = dat$Y*1000,
                            res =20e3, plot = FALSE)
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
#df$sex_complete <- unlist(lapply(strsplit(as.character(df$year_age_cohort_month_sex), " "), getElement, 5))
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

pred_out <- predict(m_age_month_cohort, newdata = pred_grid, return_tmb_object = FALSE)#, se_fit = TRUE)

# Index for entire spatial domain
years_sampled = sort(unique(hake_weight_age_df$catch_year))
pred_sims <- predict(m_age_month_cohort, newdata = pred_grid, return_tmb_object = TRUE)
index_sims = get_index(pred_sims, bias_correct = TRUE)
index_sims = index_sims %>% 
  filter(catch_year %in% years_sampled) %>% 
  mutate(region = "All")

# Compare Index between Canada and US
pred_grid_north =  pred_grid %>% 
  filter(X > 217.938)  # 217.938 is the UTM coordinate associated with 49 latitude
pred_sims_north = predict(m_age_month_cohort, newdata = pred_grid_north, return_tmb_object = TRUE)
index_sims_north = get_index(pred_sims_north)
index_sims_north = index_sims_north %>% 
  filter(catch_year %in% years_sampled) %>% 
  mutate(region = "Canada")

pred_grid_south = pred_grid %>% 
  filter(X < 217.938) 
pred_sims_south = predict(m_age_month_cohort, newdata = pred_grid_south, return_tmb_object = TRUE)
index_sims_south = get_index(pred_sims_south) 
index_sims_south = index_sims_south %>% 
  filter(catch_year %in% years_sampled) %>% 
  mutate(region = "USA")

index_sims_full = rbind(index_sims, index_sims_north, index_sims_south)

# plot the weight index by region
jpeg(filename = "plots/output_exploration/weight_index.jpeg", units="in", width=6, height=4, res = 300)
ggplot(index_sims_full, aes(x = catch_year, y = est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = region), alpha = 0.4) +
  geom_line(aes(color = region)) +
  theme_classic()
dev.off()



#save(pred_out, file = "results/sdmTMB/pred_out.RData")
pred_out %>% 
ggplot(aes(catch_year, exp(est))) + 
  geom_point() + 
  xlab("Year effects") + 
  ylab("not sure") + 
  theme_bw() 


# First run Malick's Code to create base spatial grid

year_vector = min(sort(as.numeric(unique(hake_weight_age_df$catch_year)))):max(sort(as.numeric(unique(hake_weight_age_df$catch_year))))

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

grid_pred_sdm_full2$cohort = cohort_vector[1] # add cohort
grid_pred_sdm_full3 = grid_pred_sdm_full2
for(i in 2:length(cohort_vector)) { 
  grid_pred_sdm_full2$cohort = cohort_vector[i]
  grid_pred_sdm_full3 = rbind(grid_pred_sdm_full3, grid_pred_sdm_full2)
}
