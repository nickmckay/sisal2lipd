library(tidyverse)
library(geoChronR)
library(magrittr)
library(lipdR)
fTS <- filterTs(TS,"paleoData_variable == d18O")
tts <- tidyTs(fTS,age.var = "age")

ftts <- tts %>% 
  filter(between(age,118000,130000)) %>% 
  filter(between(geo_longitude, -150,-100)) %>% 
  filter(between(geo_latitude, 0,50))

plotTimeseriesStack(ftts,time.var = "age")+theme(legend.position = "none")
mapTs(ftts)

test <- untidyTs(ftts)
pullTsVariable(test,"geo_longitude")
