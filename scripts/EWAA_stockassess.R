# load package
library(pheatmap)
library(scales)
library(colorRamps)
#source("scripts/stockassess_func.R")

#upload ewaa data from stock assessment
ewaa_obs = read.csv("raw_data/EWAA-23SA.csv", header = TRUE)
ewaa_obs = ewaa_obs[-1,] # remove year column
ewaa_obs = head(ewaa_obs, - 4) # remove last four rows with extrapolated years

# translate to wide format
ewaa_long <- pivot_longer(ewaa_obs, cols = starts_with("a"),
                   values_to = "observed", names_to = "new_age", names_prefix = "a") %>%
  arrange(as.numeric(new_age, Year)) %>% 
  rename(catch_year = Year)

# create a scaling factor for the weight values to use for alpha graphical parameter
ewaa_long <- ewaa_long %>%
  group_by(catch_year) %>% 
  mutate(N = sum(observed)) %>% 
  ungroup() %>%
  group_by(new_age) %>%
  mutate(rescale = scales::rescale(observed))
ewaa_long$new_age = as.numeric(ewaa_long$new_age)

# make heatmap of EWAA - should match stock assessment
jpeg(filename = "plots/EWAA_comparison/EWAA-heatmap.jpeg", units="in", width=8, height=7, res = 300)
ggplot(ewaa_long, aes(x = sort(as.numeric(new_age)), y = catch_year))  + 
  geom_tile(aes(alpha = rescale, fill = observed)) +
  geom_text(aes(label=round(observed,2))) +
  scale_fill_gradientn(colors = colorRampPalette(c("red",
                                                   "yellow",
                                                   "green",
                                                   "dodgerblue"))(15), guide = FALSE) +
  scale_alpha(range = c(0.1, 1)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white")) +
  labs(x = "Age", y = "Year")
dev.off()

# obtain model predictions and average across year and age to find the mean weight-at-age.add scaing factor for alpha graphical parameter. 
predicted_ewaa = combined_df %>% 
  group_by(catch_year, new_age) %>% 
  summarise(n = n(), predicted = mean(exp(est))) %>% 
  select(new_age, catch_year, predicted) %>% 
  #group_by(catch_year) %>% 
  #mutate(N = sum(predicted)) %>% 
  #ungroup() %>% 
  group_by(new_age) %>%
  mutate(rescale = scales::rescale(predicted))
  
# Make heatmap of EWAA from model predictions
jpeg(filename = "plots/EWAA_comparison/predictedEWAA-heatmap.jpeg", units="in", width=8, height=7, res = 300)
ggplot(predicted_ewaa, aes(x = as.numeric(new_age), y = catch_year))  + 
  geom_tile(aes(alpha = rescale, fill = predicted)) +
  geom_text(aes(label=round(predicted, 2))) +
  scale_fill_gradientn(colors = colorRampPalette(c("red",
                                                   "yellow",
                                                   "green",
                                                   "dodgerblue"))(15), guide = FALSE) +
  scale_alpha(range = c(0.1, 1)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"))
dev.off()


# Create new df with the deviations from the predictions and observations
ewaa_anomaly_mat = left_join(ewaa_long[,c("catch_year", "new_age","observed")], predicted_ewaa[,c("catch_year", "new_age","predicted")]) %>% 
  mutate(anom = predicted-observed)

# Make anomaly plot showing the deviations for each weight-at-age
jpeg(filename = "plots/EWAA_comparison/ewaa_anomalies.jpeg", units="in", width=8, height=7, res = 300)
ggplot(ewaa_anomaly_mat, aes(x = new_age, y = catch_year))  + 
  geom_tile(aes(fill = anom)) +
  #geom_text(aes(label=round(anom, 2))) +
  scale_fill_gradient2(na.value = "snow2") +
  #scale_alpha(range = c(0.1, 1)) +
  theme(panel.background = element_rect(fill = "white")) +
  scale_x_continuous(breaks = 0:15) +
  scale_y_continuous(breaks = 1975:2022)
dev.off()

# Look at the fishery and survey model prediction deviations independently
fishery_ewaa = combined_df %>% 
  filter(data_source == "fishery") %>% 
  group_by(catch_year, new_age) %>% 
  summarise(n = n(), predicted = mean(exp(est))) %>% 
  select(new_age, catch_year, predicted) %>% 
  #group_by(catch_year) %>% 
  #mutate(N = sum(predicted)) %>% 
  #ungroup() %>% 
  group_by(new_age) %>%
  mutate(rescale = scales::rescale(predicted))

fishery_only_anomaly = left_join(ewaa_long[,c("catch_year", "new_age","observed")], fishery_ewaa[,c("catch_year", "new_age","predicted")]) %>% 
  mutate(anom = predicted-observed)

jpeg(filename = "plots/EWAA_comparison/fishery_ewaa_anomalies.jpeg", units="in", width=8, height=7, res = 300)
ggplot(fishery_only_anomaly, aes(x = new_age, y = catch_year))  + 
  geom_tile(aes(fill = anom)) +
  #geom_text(aes(label=round(anom, 2))) +
  scale_fill_gradient2(na.value = "snow2", limits = c(-1.5,0.5)) +
  #scale_alpha(range = c(0.1, 1)) +
  theme(panel.background = element_rect(fill = "white")) +
  scale_x_continuous(breaks = 0:15) +
  scale_y_continuous(breaks = 1975:2022) +
  labs(x = "age", y = "year", title = "fishery", fill = "deviations")
dev.off()

survey_ewaa = combined_df %>% 
  filter(data_source == "survey") %>% 
  group_by(catch_year, new_age) %>% 
  summarise(n = n(), predicted = mean(exp(est))) %>% 
  select(new_age, catch_year, predicted) %>% 
  #group_by(catch_year) %>% 
  #mutate(N = sum(predicted)) %>% 
  #ungroup() %>% 
  group_by(new_age) %>%
  mutate(rescale = scales::rescale(predicted))

survey_only_anomaly = left_join(ewaa_long[,c("catch_year", "new_age","observed")], survey_ewaa[,c("catch_year", "new_age","predicted")]) %>% 
  mutate(anom = predicted-observed)

jpeg(filename = "plots/EWAA_comparison/survey_ewaa_anomalies.jpeg", units="in", width=8, height=7, res = 300)
ggplot(survey_only_anomaly, aes(x = new_age, y = catch_year))  + 
  geom_tile(aes(fill = anom)) +
  #geom_text(aes(label=round(anom, 2))) +
  scale_fill_gradient2(na.value = "snow2", limits = c(-1.5,0.5)) +
  #scale_alpha(range = c(0.1, 1)) +
  theme(panel.background = element_rect(fill = "white")) +
  scale_x_continuous(breaks = 0:15) +
  scale_y_continuous(breaks = 1975:2022) +
  labs(x = "age", y = "year", title = "survey", fill = "deviations")
dev.off()


# compare to sample size?
sample_size = read.csv("raw_data/wtatage_all_samplesize.csv")

ewaa_sample_size = sample_size %>% 
  select(-seas, -gender, -GP, -bseas, -Fleet) %>% 
  filter(Yr != -1940) %>% 
  pivot_longer(cols = starts_with("a"), values_to = "sample_size", names_to = "new_age", names_prefix = "a") %>% 
  rename(catch_year = Yr)
ewaa_sample_size$new_age = as.numeric(ewaa_sample_size$new_age)

ewaa_anomaly_mat = left_join(ewaa_anomaly_mat, ewaa_sample_size)


jpeg(filename = "plots/EWAA_comparison/deviations_by_samplesize.jpeg", units="in", width=8, height=7, res = 300)
ggplot(ewaa_anomaly_mat, aes(x = sample_size, y = anom)) +
  geom_point() +
  theme_classic() +
  geom_hline(yintercept = 0, col = "black") +
  labs(x = "sample size", y = "deviation")
dev.off()

# to make this plot, go to scripts > weight-at-age and run the first bit of lines to create the hake_df dataset
hake_df %>% 
  filter(!sex_description %in% c("Unknown/Not determined", "? Don't know where this is from. I assume this is unknown, but need to do some checking")) %>% 
  filter(catch_year %in% 2003:2015) %>% 
  mutate(points_bin = as.factor(floor(hb_latitude))) %>% 
  group_by(points_bin, sex_description, catch_year, haul_num, catch_month, catch_day, duration) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = sex_description, values_from = n) %>% 
  mutate(percent_F = Female/(Male + Female)) %>% 
  #arrange(catch_year, points_bin) %>%
  filter(!is.na(points_bin))  %>% 
  group_by(points_bin) %>% 
  mutate(sum_weights = sum(duration)) %>% 
  ungroup() %>% 
  mutate(weights = duration/sum_weights) %>% 
  group_by(points_bin) %>% 
  mutate(weights_norm = weights/mean(weights)) %>% 
  ungroup() %>% 
  ggplot(aes(x = points_bin, y = percent_F*weights_norm)) +
  geom_boxplot() +
  theme_classic() +
  geom_hline(yintercept = 0.5, lty = 2, col = "blue") +
  lims(x = as.factor(c(34:54)))
  
