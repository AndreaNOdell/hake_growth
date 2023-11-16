# Survey data ----------

bio = read.csv("raw_data/bio_info2019-2021.csv")
haul = read.csv("raw_data/haul_info2019-2021.csv")
load("results/hake_weight_age_df.RData")

# the information I need is: age, weight, catch_year, catch_month, catch_day, length, distance_fished, sex_description, hb_latitude and hb_longitude

haul_subdf = haul %>% 
  select(survey_num, haul_num, haul_back_date_time, hb_latitude, hb_longitude, haul_distance_fished_nm) %>% 
  tidyr::separate(col = haul_back_date_time , into = c('date', 'time'), sep = " ", remove = FALSE) %>% 
  tidyr::separate(col = date, into = c("catch_month", 'catch_day', "catch_year"), sep = "/") %>% 
  select(survey_num, haul_num, catch_month, catch_year, catch_day, hb_latitude, hb_longitude, haul_distance_fished_nm) %>% 
  rename(distance_fished = haul_distance_fished_nm)

bio_subdf = bio %>% 
  filter(species_name == "Pacific hake - Merluccius productus") %>% 
  select(survey_num, haul_num, sex, length_cm, weight_kg, age) %>% 
  rename(sex_description = sex, length = length_cm, weight = weight_kg) %>% 
  filter(!is.na(age)) %>% 
  mutate(new_age = age)
bio_subdf$new_age = as.integer(lapply(bio_subdf$new_age, function(x) ifelse(x > 14, 15, x)))

df = left_join(bio_subdf, haul_subdf) %>% 
  mutate(cohort = as.numeric(catch_year) - new_age) %>% 
  select(colnames(hake_weight_age_df)) %>% 
  mutate(hb_longitude = -hb_longitude)
df$catch_day = as.numeric(df$catch_day)
df$catch_month = as.numeric(df$catch_month)

hake_weight_age_df_updated = rbind(df, hake_weight_age_df)
hake_weight_age_df_updated = hake_weight_age_df_updated[complete.cases(hake_weight_age_df_updated),]

save(hake_weight_age_df_updated, file = "results/hake_weight_age_df_updated.RData")



hist(hake_weight_age_df_updated$catch_month)
plot(-hake_weight_age_df_updated$hb_longitude, hake_weight_age_df_updated$hb_latitude)


# Fishery data ------------------
us_df = read.csv("raw_data/us-weight.csv") %>% 
  select(!c(X,Length_cm)) %>% 
  rename(weight = Weight_kg, sex_description = Sex, catch_month = Month, catch_year = Year) %>% 
  mutate(new_age = ifelse(Age_yrs > 14, 15, Age_yrs))

us_df_full = us_df %>% 
  filter(sex_description %in% c("F", "M")) %>% 
  mutate(cohort = catch_year - new_age)

(nrow(us_df) - nrow(us_df_full))/nrow(us_df) # less than 2% of data unsexed/removed


can_df = read.csv("raw_data/can-weight-at-age.csv") %>% 
  select(!X) %>% 
  rename(weight = Weight_kg, sex_description = Sex, catch_month = Month, catch_year = Year) %>% 
  mutate(new_age = ifelse(Age_yrs > 14, 15, Age_yrs)) %>% 
  mutate(cohort = catch_year - new_age)

save(us_df_full, file = "results/us_df_full.RData")
save(can_df, file = "results/can_df.RData")



