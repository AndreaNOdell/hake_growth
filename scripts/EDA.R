library(tidyverse)
library(PNWColors)
library(bbmle)
library(FSA)

df = read.csv("raw_data/selection.csv") # load in data
#colnames(df)

hake_df = df %>% # make new dataframe with relevant information only
  filter(scientific_name == "Merluccius productus") %>% 
  select("age", "avg_weight", "distance_fished", #"hb_date", 
         "length", "sex_description", "weight", "catch_modified_date",
         "eq_date", "hb_date", 'hb_latitude', 'hb_longitude') %>% 
  separate(col = eq_date , into = c('date', 'time'), sep = " ", remove = FALSE) %>% 
  separate(col = date, into = c("catch_year", "catch_month", 'catch_day'), sep = "-")

# Set missing catch_year values to the year inputed for the other date related columns
hake_df[is.na(hake_df$catch_month),]$catch_year <- 2017

# Graph the observations per year
pal = pnw_palette(name="Sunset2",n=3,type="discrete") # set color palette

jpeg(file="plots/observations_per_year.jpeg")
hake_df %>% #subset data and plot observations per year
  group_by(catch_year, sex_description) %>% 
  summarise(count = n()) %>% 
  ggplot(aes(x = catch_year, y = count, fill = sex_description)) +
    geom_col() +
    theme_classic() +
    labs(title = "Number of observations per year") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_fill_manual(values = pal, limits = c('Female', 'Male', 'Unknown/Not determined'))
dev.off()

count(hake_df[hake_df$sex_description == "Female",]) # 46,816 female observations
count(hake_df[hake_df$sex_description == "Male",]) # 44,982 male observations

# Graph log-transformed weight x length data and color by sex
jpeg(file="plots/visualize_data.jpeg")
ggplot(hake_df, aes(x = log10(length), y = log10(weight), col = sex_description)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme_classic() +
  scale_colour_manual(values = pal, limits = c('Female', 'Male', 'Unknown/Not determined'))
dev.off()

####
# Length - weight model:   W = aL^b -----------------------------------------------------
####

# Fit a linear model to the log transformed data using least squares method
# log10(W) = log10(a) + b*log10(L)

# let's run the model
fit_hake_df <- hake_df[complete.cases(hake_df[ , c('weight', 'length')]), ] #subset dataset for complete cases 
growth_fit = lm(log10(weight)~log10(length),data=fit_hake_df) # run the model

# Now make predicted values
len = seq(1,100, by = 0.2) # lengths to predict
nd <- data.frame(length=len) # dataframe of lengths to predict
predicted_logW <- predict(growth_fit, nd) # predict
cf_f = exp(((log(10, base = exp(1)) * 0.04389)^2)/2) # bias correction factor (0.04389 is Residual standard error)
predicted_W = cf_f*(10^predicted_logW) # transformed back using correction factor

#plot the predictions over data and plot the residuals
jpeg(file="plots/fitted_growth_overdata.jpeg")
plot(fit_hake_df$length, fit_hake_df$weight, type = "p", main = 'Predicted growth over raw data', xlab = 'length', ylab = "weight")
lines(len, predicted_W, col = "red")
dev.off()

jpeg(file="plots/residuals.jpeg")
plot(residuals(growth_fit), type = "p", main = 'model residuals', xlab = 'residuals', ylab = "index")
abline(h = 0, col = 'red')
dev.off()

# merge the residual data back to the subsetted dataset 'fit_hake_df'
fit_hake_df = fit_hake_df %>% 
  mutate(resids = residuals(growth_fit))

# summarise growth anomalies per year
avg_resids_year = fit_hake_df %>%  #summarise residual information by year
  group_by(catch_year) %>% 
  summarise(avg = mean(resids), sd = sd(resids), n = n()) 


jpeg(file="plots/avg_residual_year.jpeg")
ggplot(avg_resids_year, aes(x = catch_year, y = avg)) +
  geom_point(aes(size = n)) +
  scale_size(range = c(1, 3)) +
  theme_classic() +
  labs(title = 'average residual growth per year', x = 'year', y = 'residual') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = 'dashed') #+
  #geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd))
dev.off()

#function to create data summary for violoin plot
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

jpeg(file="plots/avg_residual_violin_rm.jpeg")
# Violin plot of avg growth anomaly per year
ggplot(fit_hake_df, aes(x = catch_year, y = resids)) +
  geom_violin() + 
  stat_summary(fun.data=data_summary) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = 0, lty = 2) #+ 
  #coord_cartesian(ylim=c(-0.3, 0.2))
dev.off()
    # Years 1995 and 2017 have high variability in growth anomaly - potentially a response to marine heatwaves
    # I believe the variability for 2007 is due to measurement error.


fit_hake_df %>%  
  group_by(catch_year, resids) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = catch_year, y = resids, size = n)) +
    geom_point() +
    theme_classic() +
    lims(y = c(-0.5, 0.5))



# How growth anomalies vary spatially and by year?
jpeg(file="plots/spatial_resids_year.jpeg")
ggplot(fit_hake_df, aes(x = hb_longitude, y = hb_latitude, col = resids)) +
  geom_point() +
  scale_colour_gradient2(limits = c(-0.5, 0.5), mid = NA) +
  theme_classic() +
  facet_wrap(vars(catch_year)) +
  scale_x_reverse()
dev.off()

# [fit_hake_df$resids < -0.1 | fit_hake_df$resids > 0.1,]  
# make values closest to zero transparent rather than white
    # For spatial, maybe I could make small-ish bins for the latitude and longitude,
    # then summarise the information in the bin to get an average growth anomaly for the bin
  

ggplot(fit_hake_df, aes(x = hb_longitude, y = hb_latitude, col = sex_description)) +
  geom_point(alpha = 0.1) +
  #scale_colour_gradient2(limits = c(-0.5, 0.5), mid = NA) +
  theme_classic() +
  facet_wrap(vars(catch_year)) +
  scale_x_reverse() +
  scale_colour_manual(values = pal, limits = c('Female', 'Male'))


# How do growth anomalies vary by age
jpeg(file="plots/resids_age_year.jpeg")
ggplot(fit_hake_df, aes(x = age, y = resids)) +
  geom_point() +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  facet_wrap(vars(catch_year))
dev.off()

jpeg(file="plots/anomaly_timeseries_age.jpeg")
fit_hake_df %>% 
  group_by(catch_year, age) %>% 
  summarise(avg = mean(resids), sd = sd(resids), n = n()) %>% 
  ggplot(aes(x = catch_year, y = avg)) +
    geom_point() +
    theme_classic() +
    facet_wrap(vars(age)) +
    geom_hline(yintercept = 0, lty = 2) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

# the SD of the temporal growth anomalies for each age.
fit_hake_df %>% 
  group_by(age) %>% 
  summarise(avg = mean(resids), sd = sd(resids), n = n()) %>% 
  ggplot(aes(x = age, y = sd)) +
    geom_point(aes(size = n)) +
    theme_classic()



# Look at stock assessment and determiner plus age group - age 15+ in assessment
# What are the consequences of the assumptions we make about plus group

# finer temporal resolution - are fish heavier later in the season?
# fatter in the north or sampled later?

# Fit age at length models to assign ages to those that don't have them

# relative importance of the different patterns that I'm seeing - statistical model

# Is there a spatial variation in which sex dominates? 


####
# Estimating missing ages -----------------------------------------------------
####

# 2011 stock assessment fit Von Bertalanffy models 
# but they first split the data into to regimes - 1975-1989 and 1990-2010




