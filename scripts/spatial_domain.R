library(rgdal)
source("malick_code/functions.R")
load("results/hake_weight_age_df.RData")

# convert coordinates from lat long to UTM
cord.dec = SpatialPoints(cbind(-hake_weight_age_df$hb_longitude, hake_weight_age_df$hb_latitude), 
                         proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32610"))
cord.UTM <- as.data.frame(cord.UTM)

# create new columns with UTM coordinates
hake_weight_age_df$northing_m = cord.UTM[,2]
hake_weight_age_df$easting_m = cord.UTM[,1]

## Hake biomass example ------------------------------------

## 1. Create rectangular grid over entire domain
grid_pred_full <- pred_grid(x = hake_weight_age_df$easting_m, 
                            y = hake_weight_age_df$northing_m,
                            res = 30e3, plot = FALSE)
names(grid_pred_full) <- c("easting_m", "northing_m")
grid_pred_full$cell_id <- 1:nrow(grid_pred_full)
plot(grid_pred_full$easting_m, grid_pred_full$northing_m)



## 2. Exclude grid cell that > 10 km from a hake biomass location
## This can be slow b/c it needs to loop over all biomass points
excl <- exclude_grid(grid_pred_full$easting_m,
                     grid_pred_full$northing_m,
                     hake_weight_age_df$easting_m,
                     hake_weight_age_df$northing_m,
                     method = "euclidean",
                     dist = 30e3,
                     plot = FALSE)

grid_pred <- grid_pred_full[!excl, ]
row.names(grid_pred) <- NULL
save(grid_pred, file = "results/grid_pred.Rdata")


## Plot grid
plot(grid_pred$easting_m, grid_pred$northing_m)




