# Predict WAA for 1) fishery and survey without spatial information and 2) just survey with spatial information
# model includes age, year, cohort, and sex. 


# Set up ----------------
library(tidyverse)
library(sf)
library(sp)
library(rgdal)
library(ggplot2)
library(sdmTMB)
library(visreg)


# Load in Data -----------

# Survey data
load("results/hake_weight_age_df_updated.RData")
hake_sdmTMB_df = hake_weight_age_df_updated %>% 
  dplyr::select(new_age, weight, catch_year, catch_month, cohort, length, distance_fished, sex_description, hb_longitude, hb_latitude) # select relevant rows
hake_sdmTMB_df$hb_longitude = -1 * hake_sdmTMB_df$hb_longitude # fix longitude (needs to be negative)
hake_sdmTMB_df = st_as_sf(hake_sdmTMB_df, coords = c("hb_longitude", "hb_latitude")) # set spatial class
st_crs(hake_sdmTMB_df) = 4269 # assign current coordinate system (Lat/Long)
hake_sdmTMB_df = st_transform(hake_sdmTMB_df, 32610) # transform coordinate system to UTMs 
coords =  st_coordinates(hake_sdmTMB_df) 
coords = coords/1000 #convert coordinates from m to km
hake_sdmTMB_df = as.data.frame(cbind(hake_sdmTMB_df, coords, 
                                     hb_latitude = hake_weight_age_df_updated$hb_latitude,
                                     hb_longitude = -hake_weight_age_df_updated$hb_longitude)) %>%  # merge coordinates back onto dataset
  mutate(sex_description = replace(sex_description, sex_description == "Unknown/Not determined", "Unsexed"))
hake_sdmTMB_df$fcatch_year = as.factor(hake_sdmTMB_df$catch_year)
hake_sdmTMB_df$fcohort = as.factor(hake_sdmTMB_df$cohort)
hake_sdmTMB_df$sex_description = as.factor(hake_sdmTMB_df$sex_description)

survey_df = hake_sdmTMB_df

# Fishery data
load("results/can_df.RData"); load("results/us_df_full.RData") # load in cleaned data
# data for both fisheries was collected from the Hake Assessment Github > data-tables > length-weight-age > LWAdata_1975to2021.csv accessed August 2023

can_df = can_df %>% 
  mutate(country = "canada") # Canada Fishery
us_df_full = us_df_full %>% 
  mutate(country = "usa") # US Fishery
fishery_df = rbind(can_df,us_df_full) %>% # combine US and Canadian fishery
  mutate(sex_description = replace(sex_description, sex_description == "U", "Unsexed")) %>% 
  mutate(sex_description = replace(sex_description, sex_description == "M", "Male")) %>% 
  mutate(sex_description = replace(sex_description, sex_description == "F", "Female"))
fishery_df$fcatch_year = as.factor(fishery_df$catch_year)
fishery_df$sex_description <- as.factor(fishery_df$sex_description)
fishery_df$cohort <- as.integer(fishery_df$cohort)
fishery_df$fcohort <- as.factor(fishery_df$cohort)
fishery_df$country <- as.factor(fishery_df$country)


# EWAA 1  -----------------------------
# fishery and survey data, no spatial information

# merge data sets
ewaa1_df_survey = survey_df %>% # survey data
  select(new_age, weight, catch_year, fcatch_year, catch_month, cohort, sex_description, fcohort) %>% # select relevant columns
  mutate(source = "survey")
ewaa1_df_fishery = fishery_df %>% #fishery data
  select(new_age, weight, catch_year, fcatch_year, catch_month, cohort, sex_description, fcohort) %>% 
  mutate(source = "fishery")
ewaa1_df = rbind(ewaa1_df_survey, ewaa1_df_fishery) # merge data sets - 225,163 obs.
ewaa1_df$catch_year = as.numeric(ewaa1_df$catch_year)


# model
m1 = sdmTMB( 
  data = ewaa1_df,
  formula = weight ~ 1 + s(new_age) + (1|fcohort) + (1|fcatch_year) + sex_description,
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  control = sdmTMBcontrol(newton_loops = 1)
)

m1.simple = sdmTMB( 
  data = ewaa1_df,
  formula = weight ~ 0 + s(new_age),
  time_varying = ~ 1,
  time = "catch_year",
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "off",
  control = sdmTMBcontrol(newton_loops = 1)
)

# save model info
model = m1.simple
data = model$data

# expand prediction grid
pred.grid = expand.grid(new_age = 0:15,
                             catch_year = as.numeric(unique(data$catch_year))) #,
                             #sex_description = c("Male", "Female")) %>%
  #mutate(fcohort = as.factor(catch_year - new_age)) %>% 
  #mutate(fcatch_year = as.factor(catch_year))

# get predictions
preds = predict(model, newdata = pred.grid) 
preds = preds %>% 
  mutate(est_weight = exp(est))

ggplot(preds[preds$sex_description == "Female",], aes(x = new_age, y = est_weight)) + 
  geom_line(aes(col = fcohort)) +
  theme(legend.position = "none") +
  labs(subtitle = "female weight-at-age per cohort")

# data = data %>% 
#   mutate(est = preds[,ncol(preds)])

# make EWAA
ewaa1_2sources_nospatial = as.data.frame(preds %>% 
  group_by(catch_year, new_age) %>% 
  summarise(n = n(), pred_weight = mean(exp(est))))
#ewaa1_2sources_nospatial$catch_year = as.numeric(as.character(ewaa1_2sources_nospatial$fcatch_year))

ewaa1_long = ewaa1_2sources_nospatial %>% 
  select(catch_year, new_age, pred_weight)

ewaa1_wide = as.matrix(pivot_wider(ewaa1_long, names_from = new_age, values_from = pred_weight) %>% 
  relocate(`0`, .before = `1`))

write.csv(ewaa1_long, "outputs/nospatial_long_simple.csv", row.names=FALSE)
write.csv(ewaa1_wide, "outputs/nospatial_wide_simple.csv", row.names=FALSE)

# EWAA 2  -----------------------------
# survey data only with spatial information

# data is just the survey_df
ewaa2_df = survey_df[!survey_df$sex_description == "Unsexed",] %>% 
  select(new_age, weight, fcatch_year, catch_year, fcohort, catch_month, cohort, sex_description,  X, Y, hb_latitude, hb_longitude, geometry)
ewaa2_df$catch_year = as.numeric(ewaa2_df$catch_year)
ewaa2_df = droplevels(ewaa2_df)


# model
mesh <- make_mesh(ewaa2_df, xy_cols = c("X", "Y"), cutoff = 40)
  
m2 = sdmTMB( 
  data = ewaa2_df,
  formula = weight ~ 1 + s(new_age) + (1|fcohort)  + sex_description ,
  mesh = mesh,
  family = lognormal(link = "log"),
  time = "catch_year",
  spatial = "off",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016, 2018, 2020)
)

m2.simple = sdmTMB( 
  data = ewaa2_df,
  formula = weight ~ 1 + s(new_age) ,
  mesh = mesh,
  family = lognormal(link = "log"),
  time = "catch_year",
  spatial = "off",
  spatiotemporal = "ar1",
  control = sdmTMBcontrol(newton_loops = 1),
  extra_time = c(1987, 1988, 1990, 1991, 1993, 1994, 1996, 1997, 1999, 2000, 2002, 2004, 2006, 2008, 2010, 2014, 2016, 2018, 2020)
)

# save model info
model = m2.simple
data = model$data

# expand grid
pred.grid = expand.grid(new_age = 0:15, # specify all age classes
              catch_year = as.numeric(sort(unique(data$catch_year)))) #, # specify all years
              #sex_description = c("Male", "Female")) %>% # specify both sex
 # mutate(fcohort = as.factor(catch_year - new_age))
unique.coords = unique(data[,c("X", "Y")])  # , "catch_year"
  # only use spatial locations sampled in a given year (estimation crashes when estimating outside the sampling locations in a given year, but will try to resolve this later)
pred.grid.sp = cross_join(unique.coords, pred.grid)

# make predictions
preds = predict(model, newdata = pred.grid.sp)
#save(preds, file = "results/sdmTMB/ewaa_preds/preds.RData")

# quick visualization of WAA curves 
preds %>%  
  group_by(new_age, fcohort, sex_description) %>% 
  summarise(est_weight = mean(exp(est))) %>% 
  filter(sex_description == "Female") %>% 
ggplot( aes(x = new_age, y = est_weight)) + 
  geom_line(aes(col = fcohort)) +
  theme(legend.position = "none") +
  labs(subtitle = "female weight-at-age per cohort")


ewaa2_survey_spatial = as.data.frame(preds %>% 
  group_by(catch_year, new_age) %>% 
  summarise(n = n(), pred_weight = mean(exp(est))))
ewaa2_long = ewaa2_survey_spatial %>% 
  select(catch_year, new_age, pred_weight)

ewaa2_wide = as.matrix(pivot_wider(ewaa2_long, names_from = new_age, values_from = pred_weight) %>% 
  relocate(`0`, .before = `1`))

write.csv(ewaa2_long, "outputs/spatial_long_simple.csv", row.names=FALSE)
write.csv(ewaa2_wide, "outputs/spatial_wide_simple.csv", row.names=FALSE)

# EWAA 3  -----------------------------
# Fishery data only with no spatial information

# data is just the fishery_df
ewaa3_df = fishery_df %>% #fishery data
  select(new_age, weight, catch_year, fcatch_year, catch_month, cohort, sex_description, fcohort)

m3.simple = sdmTMB( 
  data = ewaa3_df,
  formula = weight ~ 0 + s(new_age),
  time_varying = ~ 1,
  family = lognormal(link = "log"),
  time = "catch_year",
  spatial = "off",
  spatiotemporal = "off",
  control = sdmTMBcontrol(newton_loops = 1)
)

# save model info
model = m3.simple
data = model$data

# expand grid
pred.grid = expand.grid(new_age = 0:15, # specify all age classes
                        catch_year = as.numeric(sort(unique(data$catch_year)))) #, # specify all years
#sex_description = c("Male", "Female")) %>% # specify both sex
# mutate(fcohort = as.factor(catch_year - new_age))

# make predictions
preds = predict(model, newdata = pred.grid)
#save(preds, file = "results/sdmTMB/ewaa_preds/preds.RData")

ewaa3_fishery = as.data.frame(preds %>% 
                                group_by(catch_year, new_age) %>% 
                                summarise(n = n(), pred_weight = mean(exp(est))))
ewaa3_long = ewaa3_fishery %>% 
  select(catch_year, new_age, pred_weight)

ewaa3_wide = as.matrix(pivot_wider(ewaa3_long, names_from = new_age, values_from = pred_weight) %>% 
                         relocate(`0`, .before = `1`))

write.csv(ewaa3_long, "outputs/fishery_only_long.csv", row.names=FALSE)
write.csv(ewaa3_wide, "outputs/fishery_only_wide.csv", row.names=FALSE)
