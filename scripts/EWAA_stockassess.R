# load package
library(pheatmap)
library(scales)
library(colorRamps)
source("scripts/stockassess_func.R")

ewaa_obs = read.csv("raw_data/EWAA-23SA.csv", header = TRUE)
head(ewaa_obs)

ewaa_obs = ewaa_obs[-1,]
ewaa_obs = head(ewaa_obs, - 4) 

#pheatmap(ewaa_obs, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, heatmap_annotations = list(alpha = lapply(seq_len(ncol(ewaa_obs)), function(col) {scale_alpha(ewaa_obs[, col], range = c(0.2, 1))})))



ewaa_long <- pivot_longer(ewaa_obs, cols = starts_with("a"),
                   values_to = "observed", names_to = "new_age", names_prefix = "a") %>%
  arrange(as.numeric(new_age, Year)) %>% 
  rename(catch_year = Year)

ewaa_long <- ewaa_long %>%
  group_by(catch_year) %>% 
  mutate(N = sum(observed)) %>% 
  ungroup() %>%
  group_by(new_age) %>%
  mutate(rescale = scales::rescale(observed))
ewaa_long$new_age = as.numeric(ewaa_long$new_age)

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

predicted_ewaa = combined_df %>% 
  group_by(catch_year, new_age) %>% 
  summarise(n = n(), predicted = mean(exp(est))) %>% 
  select(new_age, catch_year, predicted) %>% 
  #group_by(catch_year) %>% 
  #mutate(N = sum(predicted)) %>% 
  #ungroup() %>% 
  group_by(new_age) %>%
  mutate(rescale = scales::rescale(predicted))
  
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



ewaa_anomaly_mat = left_join(ewaa_long[,c("catch_year", "new_age","observed")], predicted_ewaa[,c("catch_year", "new_age","predicted")]) %>% 
  mutate(anom = predicted-observed)

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
