library("rgdal") # for `ogrInfo()` and `readOGR()`
library("tools") # for `file_path_sans_ext()`
library("dplyr") # for `inner_join()`, `filter()`, `summarise()`, and the pipe operator (%>%)
library("ggplot2") # for `fortify()` and for plotting
library("sp") # for `point.in.polygon()` and `spDists()`
library("tidyr") # for `gather()`
library("readr") # for `write_tsv()`
library("mapproj") # for ggplot:coord_map()


fortify.shape <- function(x){
  x@data$id <- rownames(x@data)
  x.f <- fortify(x, region = "id")
  x.join <- inner_join(x.f, x@data, by = "id")
}

subset.shape <- function(x, domain){
  x.subset <- filter(x, long > domain[1] & 
                       long < domain[2] & 
                       lat > domain[3] & 
                       lat < domain[4])
  x.subset
}


path.ne.coast <- ("mapping/stanford-ns372xw1938-shapefile/")
fnam.ne.coast <- "ns372xw1938.shp"
dat.coast <- readOGR(dsn = path.ne.coast, 
                     layer = file_path_sans_ext(fnam.ne.coast))
# A Large SpatialLinesDataFrame object with 4132 features and 2 fields (12.8 Mb)

# Fortify the shapefile data using `fortify.shape()`:
dat.coast <- fortify.shape(dat.coast) # a 410951x8 dataframe

# Specify the desired domain (the West Coast of the USA):
domain <- c(-135, -115, 29, 60)

# Extract the coastline data for the desired domain using `subset.shape()`:
dat.coast.wc <- subset.shape(dat.coast, domain) # a 4871x8 dataframe


# Specify the spatial extent for our map (i.e., our study area; notice that its
# dimensions are different from the domain for which we extracted the
# coastline):
xlims <- c(-135, -116)
ylims <- c(32.5, 56)

#predict_m4.st.alt = cbind(predict_m4.st.alt, 
#                          hb_longitude = hake_weight_age_df$hb_longitude, 
#                          hb_latitude = hake_weight_age_df$hb_latitude)

year_index = as.data.frame(cbind(catch_year = year_sampled, year_index = (length(year_sampled)-1):0))
hake_sdmTMB_df = left_join(hake_sdmTMB_df, year_index)
hake_sdmTMB_df = hake_sdmTMB_df %>% 
  mutate(long_for_data_map = hb_longitude - (year_index * 4))


library(PNWColors)

# Generate a base map with the coastline:
year_sampled = as.numeric(sort(unique(hake_weight_age_df_updated$catch_year)))
p0 <- ggplot() +
  #geom_hline(yintercept = 49, col = "gray") + 
  geom_point(data = pred_out[pred_out$catch_year > 1986,], aes(x = hb_longitude, y = hb_latitude, col = epsilon_st), 
             size = 1.5, shape = 15) + # , color = as.factor(sex_description)
  geom_polygon(data = dat.coast.wc, aes(x = long, y = lat, group = group), 
            color = "white", size = 0.25, fill = "grey83") + 
  coord_map(projection = "mercator") + 
  scale_x_continuous(limits = xlims, expand = c(0, 0)) + 
  scale_y_continuous(limits = ylims, expand = c(0, 0)) + 
  labs(x = "Longitude", y = "Latitude", color = "Spatiotemporal RE") +
  scale_color_gradient2(low = "darkorchid3", mid = "white",  high = "green4") +
  #scale_color_gradientn(colours = pnw_palette("Shuksan2",n = 3), values = c(-1.41,0,1)) +
  facet_wrap(~catch_year, ncol = 8) +
  theme_classic() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=12))
  

jpeg("plots/updated_output_exploration/stRE_map_year_wide.jpeg", units="in", width=15, height=9, res=300)
p0
dev.off()

jpeg("plots/output_exploration/month_effect.jpeg", units="in", width=6, height=4, res=300)
pred_out %>% 
  group_by(catch_month) %>% 
  summarise(mean_est = mean(plogis(est)), sd_est = sd(plogis(est))) %>% 
  ggplot(aes(x = catch_month, y = mean_est)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_est-sd_est, ymax = mean_est+sd_est)) +
  theme_classic() +
  labs(x = "month", y = "mean predicted weight (+/- SE)")
dev.off()

jpeg("plots/output_exploration/year_effect.jpeg", units="in", width=6, height=4, res=300)
pred_out %>% 
  group_by(catch_year) %>% 
  filter(catch_year %in% year_sampled) %>% 
  summarise(mean_est = mean(plogis(est)), sd_est = sd(plogis(est))) %>% 
  ggplot(aes(x = catch_year, y = mean_est)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_est-sd_est, ymax = mean_est+sd_est)) +
  theme_classic() +
  labs(x = "year", y = "mean predicted weight (+/- SE)")
dev.off()

jpeg("plots/output_exploration/age_effect.jpeg", units="in", width=6, height=4, res=300)
pred_out %>% 
  group_by(new_age) %>% 
  summarise(mean_est = mean(plogis(est)), sd_est = sd(plogis(est))) %>% 
  ggplot(aes(x = new_age, y = mean_est)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_est-sd_est, ymax = mean_est+sd_est)) +
  theme_classic() +
  labs(x = "age", y = "mean predicted weight (+/- SE)")
dev.off()

jpeg("plots/output_exploration/cohort_effect.jpeg", units="in", width=6, height=4, res=300)
pred_out %>% 
  group_by(cohort) %>% 
  summarise(mean_est = mean(plogis(est)), sd_est = sd(plogis(est))) %>% 
  ggplot(aes(x = cohort, y = mean_est)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_est-sd_est, ymax = mean_est+sd_est)) +
  theme_classic() +
  labs(x = "cohort", y = "mean predicted weight (+/- SE)")
dev.off()

# Q-Q plot -------
set.seed(1)
rq_res <- residuals(m_age_month_cohort) # randomized quantile residuals
jpeg("plots/output_exploration/Q-Qplot.jpeg", units="in", width=6, height=4, res=300)
qqnorm(rq_res, ylim = c(-10,10));qqline(rq_res)
dev.off()

# Weight histogram --------------
jpeg("plots/output_exploration/predicted_weight_histogram.jpeg", units="in", width=6, height=4, res=300)
hist(plogis(pred_out$est), xlab = "Weight", main = "Predicted Weights", xlim = c(0,1.5), breaks = 20)
dev.off()

jpeg("plots/output_exploration/observed_weight_histogram.jpeg", units="in", width=6, height=4, res=300)
hist(hake_weight_age_df$weight, xlim = c(0,1.5), main = "Observed Weights", breaks = 30)
dev.off()


