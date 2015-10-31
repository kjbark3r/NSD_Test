#read in data
setwd("C:/Users/kristin.barker/Documents/GitHub/NSD_Test")
data <- read.csv("Data_UTMs.csv")
data <- na.omit(data)

library(dplyr)

#remove blank locations

dailylocs <- group_by(data, Device_ID) %>%  #group all data by elk
  group_by(Date, add = TRUE) %>%  #group each elk by date
  sample_n(size = 1)  #randomly select one location (hopefully per date per elk)

write.table(dailylocs, file = "dailylocs.csv", quote = FALSE)

elklist <- distinct(data, Device_ID)
write.table(elklist, file = "elklist.csv", quote = FALSE)
