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

# NOT sex-specific
ggplot(hake_df, aes(x = log10(length), y = log10(weight), col = sex_description)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme_classic() +
  scale_colour_manual(values = pal, limits = c('Female', 'Male', 'Unknown/Not determined'))


# Length - weight model:   W = aL^b -----------------------------------------------------

# Fit a linear model to the log transformed data using least squares method
# log10(W) = log10(a) + b*log10(L)

# let's fun the model

fit_hake_df <- hake_df[complete.cases(hake_df[ , c('weight', 'length')]), ] #subset dataset for complete cases 
growth_fit = lm(log10(weight)~log10(length),data=fit_hake_df)

# Now make predicted values
len = seq(1,100, by = 0.2) # lengths to predict
nd <- data.frame(length=len) # dataframe of lengths to predict
predicted_logW <- predict(growth_fit, nd) # predict
cf_f = exp(((log(10, base = exp(1)) * 0.04389)^2)/2) # bias correction factor (0.04389 is Residual standard error)
predicted_W = cf_f*(10^predicted_logW) # transformed back

#plot the predictions over data and plot the residuals
plot(fit_hake_df$length, fit_hake_df$weight, type = "p")
lines(len, predicted_W, col = "red")
plot(residuals(growth_fit), type = "p")
abline(h = 0, col = 'red')


# merge the residual data back to the subsetted dataset 'fit_hake_df'
fit_hake_df = fit_hake_df %>% 
  mutate(resids = residuals(growth_fit))


# summarise residuals per year
avg_resids_year = fit_hake_df %>%  #summarise residual information by year and sex
  group_by(catch_year) %>% 
  summarise(avg = mean(resids), n = n()) 


ggplot(avg_resids_year, aes(x = catch_year, y = avg, size = n)) +
  geom_point() +
  scale_size(range = c(1, 3)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = 'dashed')
  





