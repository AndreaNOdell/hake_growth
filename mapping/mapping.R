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


path.ne.coast <- ("mapping/ne_10m_coastline/")
fnam.ne.coast <- "ne_10m_coastline.shp"
dat.coast <- readOGR(dsn = path.ne.coast, 
                     layer = file_path_sans_ext(fnam.ne.coast))
# A Large SpatialLinesDataFrame object with 4132 features and 2 fields (12.8 Mb)

# Fortify the shapefile data using `fortify.shape()`:
dat.coast <- fortify.shape(dat.coast) # a 410951x8 dataframe

# Specify the desired domain (the West Coast of the USA):
domain <- c(-128, -115, 29, 60)

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

# Generate a base map with the coastline:
p0 <- ggplot() + 
  geom_point(data = df[df$catch_year != 1986,], aes(x = (hb_longitude*-1), y = hb_latitude, color = omega_s), 
             size = 0.8, shape = 15) +
  geom_path(data = dat.coast.wc, aes(x = long, y = lat, group = group), 
            color = "black", size = 0.25) + 
  coord_map(projection = "mercator") + 
  scale_x_continuous(limits = xlims, expand = c(0, 0)) + 
  scale_y_continuous(limits = ylims, expand = c(0, 0)) + 
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_gradient2(low = "red", mid = "white",  high = "blue") +
  #facet_wrap(~catch_year) +
  theme_classic()
#jpeg(filename = "plots/sdmTMB/mapping/m4.st.stRE.jpeg")

jpeg("plots/nested_models/spatialREmap.jpeg", units="in", width=5, height=5, res=300)
p0
dev.off()







