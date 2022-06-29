library(tidyverse)
library(PNWColors)
library(bbmle)

df = read.csv("raw_data/selection.csv")
colnames(df)

hake_df = df %>% 
  filter(scientific_name == "Merluccius productus") %>% 
  select("age", "avg_weight", "distance_fished", #"hb_date", 
         "length", "sex_description", "weight", "catch_modified_date",
         "eq_date", "hb_date") %>% 
  separate(col = eq_date , into = c('date', 'time'), sep = " ", remove = FALSE) %>% 
  separate(col = date, into = c("catch_year", "catch_month", 'catch_day'), sep = "-")

# Set missing catch_year values to the year inputed for the other date related columns
hake_df[is.na(hake_df$catch_month),]$catch_year <- 2017

# Graph of the observations per year
pal = pnw_palette(name="Starfish",n=3,type="discrete")

hake_df %>% 
  group_by(catch_year, sex_description) %>% 
  summarise(count = n()) %>% 
  ggplot(aes(x = catch_year, y = count, fill = sex_description)) +
    geom_col() +
    theme_classic() +
    labs(main = "Number of observations per year") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_fill_manual(values = pal, limits = c('Female', 'Male', 'Unknown/Not determined'))

count(hake_df[hake_df$sex_description == "Female",]) # 46,816 female observations
count(hake_df[hake_df$sex_description == "Male",]) # 44,982 male observations

hake_df %>% 
  filter(sex_description == "Female") %>% 
  ggplot(aes(x = log(length), y = log(weight))) +
    geom_point() +
    #lims(x = c(0,85), y = c(0,4.2)) +
    #facet_wrap(vars(catch_year)) +
    theme_classic() +
    labs(title = "Hake growth - Female") #+
    #scale_colour_continuous(limits = c(0,4.2))

hake_df %>% 
  filter(sex_description == "Male") %>% 
  ggplot(aes(x = log(length), y = log(weight))) +
  geom_point() +
    #lims(x = c(0,85), y = c(0,4.2)) +
    theme_classic() +
    labs(title = "Hake growth - Male") #+
    #scale_colour_continuous(limits = c(0,4.2))

#
# Length - weight model:   W = aL^b -----------------------------------------------------

# Fit a linear model to the log transformed data using least squares method
# log10(W) = log10(a) + b*log10(L)

# let's fun the model
female_hake_df = hake_df %>% 
  filter(sex_description == "Female") 
fit_female_df <- female_hake_df[complete.cases(female_hake_df[ , c('weight', 'length')]), ] #subset dataset for complete cases 
female_fit = lm(log10(weight)~log10(length),data=fit_female_df)

male_hake_df = hake_df %>% 
  filter(sex_description == "Male")
fit_male_df <- male_hake_df[complete.cases(male_hake_df[ , c('weight', 'length')]), ] #subset dataset for complete cases 
male_fit = lm(log10(weight)~log10(length),data=fit_male_df)

# Now make predicted values
len = seq(1,100, by = 0.2)
nd <- data.frame(length=len)
#females
predicted_logW_F <- predict(female_fit, nd)
cf_f = exp(((log(10, base = exp(1)) * 0.04327)^2)/2)
predicted_W_F = cf_f*(10^predicted_logW_F) # transformed back
#males
predicted_logW_M <- predict(male_fit, nd)
cf_m = exp(((log(10, base = exp(1)) * 0.04183)^2)/2)
predicted_W_M = cf_m*(10^predicted_logW_M)

#plot the predicted against observed
plot(female_hake_df$length, female_hake_df$weight, type = "p")
lines(len, predicted_W_F, col = "red")
plot(residuals(female_fit), type = "p")
abline(h = 0, col = 'red')
plot(male_hake_df$length, male_hake_df$weight, type = "p")
lines(len, predicted_W_M, col = "red")


# let's merge the residual data back to the subsetted dataset 'fit_sex_df'
fit_female_df = fit_female_df %>% 
  mutate(resids = residuals(female_fit))

fit_male_df = fit_male_df %>% 
  mutate(resids = residuals(male_fit))


# summarise residuals per year
avg_resids_year_f = fit_female_df %>%  #summarise residual information by year and sex
  group_by(catch_year) %>% 
  summarise(avg = mean(resids), n = n()) %>% 
  mutate(sex = "female")

avg_resids_year_m = fit_male_df %>%  #summarise residual information by year and sex
  group_by(catch_year) %>% 
  summarise(avg = mean(resids), n = n()) %>% 
  mutate(sex = "male")

ave_resids_year = rbind(avg_resids_year_f, avg_resids_year_m)  #Combine male and female

ggplot(ave_resids_year, aes(x = catch_year, y = avg, col = sex, size = n)) +
  geom_point() +
  scale_size(range = c(1, 3)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 'dashed')
  





