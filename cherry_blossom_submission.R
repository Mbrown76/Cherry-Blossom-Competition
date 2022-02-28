###########################################################
  # Marek Brown
  # Trenton Kealey
  # Kirsten Rice

  # Cherry Blossom Prediction Code

  # February 28, 2022
###########################################################

# Prep and load packages
library(tidyverse)
library(readr)
library(Metrics)
library(mgcv)
library(rio)

# load in peak bloom date data
cherry <- read.csv("data/washingtondc.csv") %>% 
  bind_rows(read.csv("data/liestal.csv")) %>% 
  bind_rows(read.csv("data/kyoto.csv"))

# load in temperature data
dctemp <- read.csv("data/dcatemps.csv")
liestaltemp <- read.csv("data/Liestal_Temps_Spring.csv")
kyototemp <- read.csv("data/kyoto_temps.csv")

# dataset formatting

dctemp <- dctemp[,c(1:14,16)] # gets rid of extra dc variables
dctemp$location <-"washingtondc" # adds location to dataset

colnames(dctemp) <- c("year", "Jan", "Feb", 'Mar', 'Apr', 'May', 'Jun', 'Jul',
                      'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Annual', 'SPRING', 'location') 
# converting dc temperatures from Fahrenheit to Celcius to match Kyoto and Liestal
dctemp <-
  dctemp %>%
  mutate(
    Jan = (5/9) * (Jan - 32), Jan = round(Jan, digits = 2),
    Feb = (5/9) * (Feb - 32), Feb = round(Feb, digits = 2),
    Mar = (5/9) * (Mar - 32), Mar = round(Mar, digits = 2),
    Apr = (5/9) * (Apr - 32), Apr = round(Apr, digits = 2),
    May = (5/9) * (May - 32), May = round(May, digits = 2),
    Jun = (5/9) * (Jun - 32), Jun = round(Jun, digits = 2),
    Jul = (5/9) * (Jul - 32), Jul = round(Jul, digits = 2),
    Aug = (5/9) * (Aug - 32), Aug = round(Aug, digits = 2),
    Sep = (5/9) * (Sep - 32), Sep = round(Sep, digits = 2),
    Oct = (5/9) * (Oct - 32), Oct = round(Oct, digits = 2),
    Nov = (5/9) * (Nov - 32), Nov = round(Nov, digits = 2),
    Dec = (5/9) * (Dec - 32), Dec = round(Dec, digits = 2),
    SPRING = (5/9) * (SPRING - 32), SPRING = round(SPRING, digits = 2),
    Annual = (5/9) * (Annual - 32), Annual = round(Annual, digits = 2),
  )

# fix one date and add location for Kyoto
kyototemp[1,1] <- 1881
for (i in 1:141){
  kyototemp$location[i] = "kyoto"
}

# Change column names in each temperature dataset to be consistent
colnames(kyototemp) <- c("year", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", 'Aug',
                         "Sep", "Oct", "Nov", "Dec", "Annual", "SPRING", "location")



colnames(liestaltemp) <- c("year", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", 'Aug',
                           "Sep", "Oct", "Nov", "Dec","SPRING", "Annual", "location")



# combine temp datasets
dcjapanswitz <-
  rbind(dctemp,kyototemp,liestaltemp)
view(dcjapanswitz)

# join temperature and bloom doy datasets
dsj_cherry <-
  dcjapanswitz %>%
  left_join(cherry, by = c("year", "location" ))

# gets rid of years where there is temperature data but no bloom doy data
dsjcherry <- na.omit(dsj_cherry)

# model building
annreg <- lm(bloom_doy~Annual+long+lat, data=dsjcherry)
summary(annreg)
par(mfrow=c(2,2))
plot(annreg)

fm_poly <- lm(bloom_doy~Feb+poly(Mar,4)+lat+long, data=dsjcherry)
summary(fm_poly)
par(mfrow=c(2,2))
plot(fm_poly)

### model validation
# validation by mae and rmse
m <- mae(dsjcherry$bloom_doy, fm_poly$fitted.values)
r <- rmse(dsjcherry$bloom_doy, fm_poly$fitted.values)
m
r
r/m
m2 <- mae(dsjcherry$bloom_doy, annreg$fitted.values)
r2 <- rmse(dsjcherry$bloom_doy, annreg$fitted.values)
m2
r2
r2/m2

# validation by LOOCV
oos_preds <- numeric(nrow(dsjcherry))
for(i in 1:nrow(dsjcherry)) {
  model_sub <- update(object = fm_poly, data = dsjcherry[-i,])
  oos_preds[i] <- abs(dsjcherry$bloom_doy[i] - predict(model_sub, newdata = dsjcherry[i, ]))
}

mean(oos_preds)


# validation by location
model_sub <- update(object = fm_poly, data = subset(dsjcherry, location %in% c('kyoto', 'liestal')))
oos_preds_dc <- abs(dsjcherry$bloom_doy - predict(model_sub, newdata = subset(dsjcherry, location = 'washingtondc')))
mean(oos_preds_dc)

model_sub <- update(object = fm_poly, data = subset(dsjcherry, location %in% c('kyoto', 'washingtondc')))
oos_preds_l <- abs(dsjcherry$bloom_doy - predict(model_sub, newdata = subset(dsjcherry, location = 'liestal')))
mean(oos_preds_l)

model_sub <- update(object = fm_poly, data = subset(dsjcherry, location %in% c('washingtondc', 'liestal')))
oos_preds_k <- abs(dsjcherry$bloom_doy - predict(model_sub, newdata = subset(dsjcherry, location = 'kyoto')))
mean(oos_preds_k)

oos_preds_all <- c(oos_preds_dc,oos_preds_l,oos_preds_k)
mean(oos_preds_all)


# Predictions

# DC PREDICTIONS

fut_years <- data.frame(year=c(2023:2032)) # set data for 10 years after 2022
dc_data <-dsjcherry %>%
  filter(location == 'washingtondc', year > 1961) # we want to predict with last 50 years only

dc_temp_reg <- gam(Annual~year, data=dc_data) # local regression to predict future annual temperature
summary(dc_temp_reg)
dc_temp_pred <- predict(dc_temp_reg, newdata=fut_years) # use regression to predict temperatures
dc_temp_pred
dc_temp_predictions <- expand_grid(year = 2023:2032) %>% 
  bind_cols(predicted_annual_temp = dc_temp_pred) # bind predictions with year

dc_lat <- c(38.88535, 38.88535, 38.88535, 38.88535,38.88535,
            38.88535,38.88535,38.88535,38.88535, 38.88535)
dc_long <- c(-77.03863, -77.03863, -77.03863, -77.03863,
             -77.03863, -77.03863, -77.03863, -77.03863, 
             -77.03863, -77.03863)
dc_latlong <- cbind(dc_lat, dc_long) # set latitude and longitude for 2022 predictions

dc_temp_predictions <- cbind(dc_temp_predictions,dc_latlong) # bind predictions with location
colnames(dc_temp_predictions) <- c("year", "Annual", "lat", "long") # rename columns
dc_pred2 <- predict(annreg, newdata=dc_temp_predictions) # make predictions for bloom doy using predicted temperature
dc_final <- cbind(dc_temp_predictions[,1], dc_pred2) # bind temperature with predictions
colnames(dc_final) <- c("year", "predicted_doy")

dc_2022 <- data.frame(5, 9.44, 38.88535, -77.03863) # set 2022 variables
colnames(dc_2022) <- c('Feb', 'Mar', 'lat', 'long') # format names
dc_2022_pred <- data.frame(2022, predict(fm_poly, newdata=dc_2022)) # predict on 2022 variable values
colnames(dc_2022_pred) <- c('year', 'predicted_doy') # format names
dc_final <- rbind(dc_2022_pred, dc_final) # bind 2022 predictions with 2023-2032 predictions
dc_final <- # show and bind year, bloom doy, and that year's date
  dc_final %>% 
  mutate(predicted_doy = round(predicted_doy),
         predicted_date = strftime(
           strptime(sprintf('%04d-01-01', year), '%Y-%m-%d') + 
             as.difftime(predicted_doy, units = 'days'),
           '%Y-%m-%d'))
# 2022 temperatures were found and calculated outside of R. More information at bottom of document
# All other locations follow this exact method with different variable names and values


# LIESTAL PREDICTIONS

fut_years <- data.frame(year=c(2023:2032))
liestal_data <-dsjcherry %>%
  filter(location == 'liestal', year > 1961)

l_temp_reg <- gam(Annual~year, data=liestal_data) 
summary(l_temp_reg)
l_temp_pred <- predict(l_temp_reg, newdata=fut_years)
l_temp_pred
l_temp_predictions <- expand_grid(year = 2023:2032) %>% 
  bind_cols(predicted_annual_temp = l_temp_pred)

liestal_lat <- c(47.48140, 47.48140, 47.48140, 47.48140, 47.48140,
                 47.48140, 47.48140, 47.48140, 47.48140, 47.48140)
liestal_long <- c(7.730519, 7.730519, 7.730519, 7.730519, 7.730519,
                  7.730519, 7.730519, 7.730519, 7.730519, 7.730519)
liestal_latlong <- cbind(liestal_lat, liestal_long)

l_temp_predictions <- cbind(l_temp_predictions,liestal_latlong)
colnames(l_temp_predictions) <- c("year", "Annual", "lat", "long")
l_pred2 <- predict(annreg, newdata=l_temp_predictions)
liestal_final <- cbind(l_temp_predictions[,1], l_pred2)
colnames(liestal_final) <- c("year", "predicted_doy")

l_2022 <- data.frame(5.5, 4.45, 47.48140, 7.730519)
colnames(l_2022) <- c('Feb', 'Mar', 'lat', 'long')
l_2022_pred <- data.frame(2022, predict(fm_poly, newdata=l_2022))
colnames(l_2022_pred) <- c('year', 'predicted_doy')
liestal_final <- rbind(l_2022_pred, liestal_final)
liestal_final <- 
  as.data.frame(liestal_final) %>% 
  mutate(predicted_doy = round(predicted_doy),
         predicted_date = strftime(
           strptime(sprintf('%04d-01-01', year), '%Y-%m-%d') + 
             as.difftime(predicted_doy, units = 'days'),
           '%Y-%m-%d'))
# 2022 temperatures were found and calculated outside of R. More information at bottom of document
# All other locations follow this exact method with different variable names and values

# KYOTO PREDICTIONS

fut_years <- data.frame(year=c(2023:2032))
kyoto_data <-dsjcherry %>%
  filter(location == 'kyoto', year > 1961)

k_temp_reg <- gam(Annual~year, data=kyoto_data) 
summary(k_temp_reg)
k_temp_pred <- predict(k_temp_reg, newdata=fut_years)
k_temp_pred
k_temp_predictions <- expand_grid(year = 2023:2032) %>% 
  bind_cols(predicted_annual_temp = k_temp_pred)

kyoto_lat <- c(35.01198, 35.01198, 35.01198, 35.01198, 35.01198,
               35.01198, 35.01198, 35.01198, 35.01198, 35.01198)
kyoto_long <- c(135.67611, 135.67611, 135.67611, 135.67611, 135.67611,
                135.67611, 135.67611, 135.67611, 135.67611, 135.67611)
kyoto_latlong <- cbind(kyoto_lat, kyoto_long)

k_temp_predictions <- cbind(k_temp_predictions,kyoto_latlong)
colnames(k_temp_predictions) <- c("year", "Annual", "lat", "long")
k_pred2 <- predict(annreg, newdata=k_temp_predictions)
kyoto_final <- cbind(k_temp_predictions[,1], k_pred2)
colnames(kyoto_final) <- c("year", "predicted_doy")

k_2022 <- data.frame(5.06, 8.34, 47.48140, 7.730519)
colnames(k_2022) <- c('Feb', 'Mar', 'lat', 'long')
k_2022_pred <- data.frame(2022, predict(fm_poly, newdata=k_2022))
colnames(k_2022_pred) <- c('year', 'predicted_doy')
kyoto_final <- rbind(k_2022_pred, kyoto_final)
kyoto_final <- 
  as.data.frame(kyoto_final) %>% 
  mutate(predicted_doy = round(predicted_doy),
         predicted_date = strftime(
           strptime(sprintf('%04d-01-01', year), '%Y-%m-%d') + 
             as.difftime(predicted_doy, units = 'days'),
           '%Y-%m-%d'))
# 2022 temperatures were found and calculated outside of R. More information at bottom of document
# All other locations follow this exact method with different variable names and values

# VANCOUVER PREDICTIONS

fut_years <- data.frame(year=c(2023:2032))
van_data <-dsjcherry %>%
  filter(year > 1961)

van_temp_reg <- lm(Annual~year+lat+long, data=van_data) # this regression includes location
    # because all three other locations are being used to predict. this gives us a better prediction model
summary(van_temp_reg)
van_lat <- c(49.2827, 49.2827, 49.2827, 49.2827, 49.2827,
             49.2827, 49.2827, 49.2827, 49.2827, 49.2827)
van_long <- c(123.1207, 123.1207, 123.1207, 123.1207, 123.1207,
              123.1207, 123.1207, 123.1207, 123.1207, 123.1207)
van_years <- cbind(fut_years, van_lat, van_long)
colnames(van_years) <- c('year', 'lat', 'long')
van_temp_pred <- predict(van_temp_reg, newdata=van_years)

van_temp_pred
van_temp_predictions <- expand_grid(year = 2023:2032) %>% 
  bind_cols(predicted_annual_temp = van_temp_pred)

van_latlong <- cbind(van_lat, van_long)

van_temp_predictions <- cbind(van_temp_predictions,van_latlong)
colnames(van_temp_predictions) <- c("year", "Annual", "lat", "long")
van_pred2 <- predict(annreg, newdata=van_temp_predictions)
van_final <- cbind(van_temp_predictions[,1], van_pred2)
colnames(van_final) <- c("year", "predicted_doy")

van_2022 <- data.frame(5.02, 6.02, 49.2827, 123.1207)
colnames(van_2022) <- c('Feb', 'Mar', 'lat', 'long')
van_2022_pred <- data.frame(2022, predict(fm_poly, newdata=van_2022))
colnames(van_2022_pred) <- c('year', 'predicted_doy')
van_final <- rbind(van_2022_pred, van_final)
van_final <- 
  as.data.frame(van_final) %>% 
  mutate(predicted_doy = round(predicted_doy),
         predicted_date = strftime(
           strptime(sprintf('%04d-01-01', year), '%Y-%m-%d') + 
             as.difftime(predicted_doy, units = 'days'),
           '%Y-%m-%d'))

# 2022 temperatures were found and calculated outside of R. More information at bottom of document
# All other locations follow this exact method with different variable names and values



dc_final
liestal_final
kyoto_final
van_final

#################################################################################################

# Appendix - links to 2022 data information

# https://www.accuweather.com/en/us/washington/20006/february-weather/327659?year=2022
# https://www.accuweather.com/en/us/washington/20006/march-weather/327659?year=2022
# https://www.accuweather.com/en/ch/liestal/311994/february-weather/311994
# https://www.accuweather.com/en/ch/liestal/311994/march-weather/311994?year=2022
# https://www.accuweather.com/en/jp/kyoto-shi/224436/february-weather/224436
# https://www.accuweather.com/en/jp/kyoto-shi/224436/march-weather/224436?year=2022
# https://www.accuweather.com/en/ca/vancouver/v6c/february-weather/53286
# https://www.accuweather.com/en/ca/vancouver/v6c/march-weather/53286?year=2022
# ALL WEATHER LINKS ACCESSED AND CALCULATED FEB 25 2022







