####################### Calculate Degree Days and Dalight hours ##################
library(suncalc)
library(dplyr)

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data5_260309.RDATA')

#calculate Daylight hours ...............................................
plankton_RUE_data6 <- plankton_RUE_data5
sun_data <- plankton_RUE_data6[, c("Date", "latitude", "longitude")]
names(sun_data) <- c("date", "lat", "lon")
sun <- getSunlightTimes(
  data = sun_data,
  keep = c("sunrise", "sunset")
)

sun$daylight_hours <- as.numeric(
  difftime(sun$sunset, sun$sunrise, units = "hours")
)

plankton_RUE_data6$daylight_hours <- sun$daylight_hours

save(plankton_RUE_data6, file = 'plankton_RUE_data6_260309.RDATA')

