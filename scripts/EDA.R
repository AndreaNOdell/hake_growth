library(tidyverse)
library(PNWColors)
library(bbmle)
library(fishmethods)

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

jpeg(file="plots/length_freq_year.jpeg")
hake_df %>% 
  group_by(catch_year, length) %>% 
  summarise(n = n())  %>% 
  ggplot(aes(x = length, y = n)) +
  geom_col() +
  theme_classic() +
  facet_wrap(vars(catch_year))
dev.off()


# spatial variation from year to year in sampling
jpeg(file="plots/spatial_sampling_year.jpeg")
ggplot(hake_df, aes(x = hb_longitude, y = hb_latitude)) +
  geom_point() +
  #scale_colour_gradient2(limits = c(-0.5, 0.5), mid = NA) +
  theme_classic() +
  facet_wrap(vars(catch_year)) +
  scale_x_reverse() +
  labs(main = "spatial variation in sampling per year")
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
jpeg(file="plots/variable_resids_age.jpeg")
fit_hake_df %>% 
  group_by(age) %>% 
  summarise(avg = mean(resids), sd = sd(resids), n = n()) %>% 
  ggplot(aes(x = age, y = sd)) +
    geom_point(aes(size = n)) +
    theme_classic()
dev.off()



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

# I'm going to start by fitting a curve to all the data.
age_hake_df = hake_df[complete.cases(hake_df[ , c('age', 'length')]), ] # subset the data for only complete cases of age and length
growth_aggregated = growth(intype = 1, unit = 1, size = age_hake_df$length, age = age_hake_df$age, 
                      Sinf = 80, K = 0.5, t0 = 4)

Sinf_aggregated = 51.6616
K_aggregated = 0.3655
t0_aggregated = -0.7004
# RSS = 551852

jpeg(file="plots/aggregate_growth_curve.jpeg")
plot(age_hake_df$age, age_hake_df$length, type = "p", col = 'blue')
lines(0:25, Sinf_aggregated * (1 - exp(-(K_aggregated * (0:25 - t0_aggregated)))), type = 'l', lwd = 2, col = "red")
dev.off()


# Two regimes (1980-2003; 2004-2017)
pre_2003_df = age_hake_df %>%  # 24,403 observations
  filter(catch_year < 2004)

growth_pre2003 = growth(intype = 1, unit = 1, size = pre_2003_df$length, age = pre_2003_df$age, 
                           Sinf = 80, K = 0.5, t0 = 4)

Sinf_pre2003 = 51.9182
K_pre2003 = 0.3634
t0_pre2003 = -0.891
#RSS = 324002

plot(age_hake_df$age, age_hake_df$length, type = "p", col = 'blue')
lines(0:25, Sinf_pre2003 * (1 - exp(-(K_pre2003 * (0:25 - t0_pre2003)))), type = 'l', lwd = 2)

post_2003_df = age_hake_df %>%  # 24,403 observations
  filter(catch_year > 2003)

growth_post2003 = growth(intype = 1, unit = 1, size = post_2003_df$length, age = post_2003_df$age,
                        Sinf = 80, K = 0.5, t0 = 4)

Sinf_post2003 = 51.1816
K_post2003 = 0.3501
t0_post2003 = -0.7535     
# RSS = 200280

plot(age_hake_df$age, age_hake_df$length, type = "p", col = 'blue')
lines(0:25, Sinf_post2003 * (1 - exp(-(K_post2003 * (0:25 - t0_post2003)))), type = 'l', lwd = 2)


# These don't look good.....

# Maybe I should concatenate the plus group first, so that there are more observations 
# for those age classes.

# assign all ages 15+ the age 15. 
hake_df = hake_df %>% 
  mutate(new_age = age)
hake_df$new_age = as.integer(lapply(hake_df$new_age, function(x) ifelse(x > 14, 15, x)))
# fit model to complete weight and length cases
fit_hake_df <- hake_df[complete.cases(hake_df[ , c('weight', 'length')]), ]#subset dataset for complete cases 
fit_hake_df$cohort = as.numeric(fit_hake_df$catch_year) - fit_hake_df$new_age

growth_fit = lm(log10(weight)~log10(length),data=fit_hake_df) # run the model
# add residual info back to dataset
fit_hake_df = fit_hake_df %>% 
  mutate(resids = residuals(growth_fit))

# Let's take a quick look at the variability of growth anomalies by age
# i.e.
jpeg(file="plots/resid_sd_age_plusgroup.jpeg")
fit_hake_df %>% 
  group_by(new_age) %>% 
  summarise(avg = mean(resids), sd = sd(resids), n = n()) %>% 
  mutate(positive = avg > 0) %>% 
  ggplot(aes(x = new_age, y = sd, col = positive)) +
  geom_point(aes(size = n)) +
  theme_classic() +
  labs(y = "residual standard deviation", x = "age", title = "growth anomaly variability by age")
dev.off()


age_hake_df = hake_df[complete.cases(hake_df[ , c('new_age', 'length')]), ] # subset the data for only complete cases of age and length

growth_plusgroup = growth(intype = 1, unit = 1, size = age_hake_df$length, age = age_hake_df$new_age, 
                           Sinf = 80, K = 0.5, t0 = 4)

Sinf_plus = 51.6614
K_plus = 0.3655
t0_plus = -0.7002     
# RSS = 551912

plot(age_hake_df$age, age_hake_df$length, type = "p", col = 'blue')
lines(0:25, Sinf_plus * (1 - exp(-(K_plus * (0:25 - t0_plus)))), type = 'l', lwd = 2, lty = 2)


# It doesn't matter whether you keep all ages or create a plus group - the growth
# model remains basically the same

ggplot(age_hake_df, aes(x = age, y = length)) +
  geom_point() +
  theme_classic() +
  facet_wrap(vars(catch_year))


# Frequency distribution of lengths per age
jpeg(file="plots/length_freq_age.jpeg")
ggplot(hake_df, aes(x = length, col = new_age)) +
  geom_freqpoly() +
  theme_classic()
dev.off()


ggplot(hake_df, aes(x = length, col = new_age)) +
  geom_freqpoly() +
  theme_classic() +
  facet_wrap(vars(catch_year))



# Add estimated ages----------------------------------------------------------------------

plot(hake_df$age, hake_df$length, type = "p")
lines(seq(0,25,by = 0.1), 49.71 - (49.71 - 1)*exp(-0.51*seq(0,25,by = 0.1)), col = 'red')

plot(hake_df$length, hake_df$age, type = "p")
lines(x = seq(0.1,49,by = 0.1), t, col = "red")

calc_age = function(length) {
  if(length < 49.71) {
    x = (log(-c(length) + 49.71) - log(49.71 - 1)) /(-0.51)
    return(x)
  } else {
    x = 15
    return(x)
  }
}

hake_no_age <- hake_df[is.na(hake_df$age),] 
hake_no_age = hake_no_age[complete.cases(hake_no_age[ , c('length')]), ] %>% 
  mutate(new_age = if_else(length < 49.71, ((log(-length + 49.71) - log(49.71 - 1)) /(-0.51)), rlnorm(1, 2.127496758, 0.283937767))) 
hake_no_age$new_age = round(hake_no_age$new_age)

unique(hake_no_age$new_age)


# Let's look at the distribution of fish ages for fish greater than 49.71cm
library(MASS)
x = na.omit(hake_df[hake_df$length > 49.71,]$age)
fitdistr(x, "lognormal")
hist(x, prob = TRUE)
curve(dlnorm(x, 2.127496758, 0.283937767), col = 2, add = TRUE)

rlnorm(1, 2.127496758, 0.283937767)



# let's look at finer temporal resolution   ---------------------------------
head(fit_hake_df)

# this plot shows the avg resids per month per year
jpeg(file="plots/resids_month_year.jpeg")
fit_hake_df %>% 
  group_by(catch_year, catch_month) %>% 
  summarise(avg = mean(resids), sd = sd(resids), n = n()) %>% 
  filter(!is.na(catch_month)) %>% 
  ggplot(aes(x = catch_month, y = avg)) +
      geom_point(aes(size = n)) +
      geom_hline(yintercept = 0, lty = 2) +
      facet_wrap(vars(catch_year)) +
      theme_classic() +
      labs(title = "variability in growth anomalies within sampling year", x = "sampling month", y = "avg growth anomaly")
dev.off()


# this plot shows the sampling movement up the coast through the summer
fit_hake_df %>% 
  group_by(catch_year, catch_month) %>% 
  summarise(avglat = mean(hb_latitude), avgweight = mean(weight)) %>% 
  filter(!is.na(catch_month)) %>% 
  ggplot(aes(x = catch_month, y = avglat, size = avgweight)) +
      geom_point() +
      #geom_hline(yintercept = 0, lty = 2) +
      facet_wrap(vars(catch_year)) +
      theme_classic()

pal2 = pnw_palette(name="Sunset2",n=5,type="discrete") # set color palette

jpeg(file="plots/growth_variability_3.7.15.jpeg")
fit_hake_df %>% 
  group_by(catch_year, age) %>% 
  summarise(avgresids = mean(resids)) %>% 
  filter(age %in% c(3,7,15)) %>% 
  ggplot(aes(x = as.factor(catch_year), y = avgresids, colour = as.factor(age), group = as.factor(age))) +
      geom_line() +
      geom_hline(yintercept = 0, lty = 2) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      labs(title = "Variability in growth anomalies through time", subtitle = "for ages 3, 7, and 15", x = "Year", y = "avg growth anomaly") +
      scale_colour_manual(values = pal2)
dev.off()


# Looking at average weight per month per year
jpeg(file="plots/weight_month_year.jpeg")
fit_hake_df %>% 
  group_by(catch_year, catch_month) %>% 
  summarise(avg = mean(weight), sd = sd(resids), n = n()) %>% 
  filter(!is.na(catch_month)) %>% 
  ggplot(aes(x = catch_month, y = avg)) +
  geom_point(aes(size = n)) +
  #geom_hline(yintercept = 0, lty = 2) +
  facet_wrap(vars(catch_year)) +
  theme_classic() +
  labs(title = "variability in weight within sampling year", x = "sampling month", y = "avg weight")
dev.off()


jpeg(file="plots/condition_variability_cohort.jpeg")
fit_hake_df %>% 
  group_by(catch_year, cohort) %>% 
  summarise(avgresids = mean(resids), n = n()) %>% 
  filter(cohort %in% c(1993, 1996, 1999, 2005, 2010, 2014)) %>% 
  ggplot(aes(x = as.factor(catch_year), y = avgresids, colour = as.factor(cohort), group = as.factor(cohort))) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title = "Variability in weight-length anomalies through time", subtitle = "for cohort", x = "Year", y = "avg growth anomaly") #+
  #scale_colour_manual(values = pal2)
dev.off()




