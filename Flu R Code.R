setwd("C:/Users/icamo/OneDrive/Documents/Clemson/Fall 2018/MATH 8820/Project 1")
library(coda)
library(tidyr)
library(MASS)

precip <- read.csv("Precipitation.csv", header= T, fileEncoding="UTF-8-BOM")
precip <- precip[precip$State != "FL", ]
precip <- precip[,1:14]


states <- droplevels(unique(precip$State))

temp <- read.csv("Temperature.csv", header= T, fileEncoding="UTF-8-BOM")
temp <- temp[temp$State != "FL", ]

unemp <- read.csv("Imputed_Unemployment_Data.csv", header = T, fileEncoding = "UTF-8-BOM")
colnames(unemp) <- c('Year', 'Month', as.character(states))

pop.land <- read.csv("PopulationLandArea.csv", header= T, fileEncoding="UTF-8-BOM")

pop.land <- pop.land[,c(5,7,8,9,10,11,12,13,14,15,58)]
pop.land <- pop.land[c(-2,-11),]
pop.land <- pop.land[pop.land$NAME != "Florida", ]

land.area <- pop.land[,11]
pop.density <- pop.land[,2:10]/land.area
pop.density <- cbind(pop.land$NAME, pop.density)

colnames(pop.density) <- c("State", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018")

temp <- temp %>%
  gather('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', key = "Month", value = "Deg.Far")
map = setNames(1:12, c('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'))
temp[,3] <- map[unlist(temp[,3])]

precip <- precip %>%
  gather('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', key = "Month", value = "Inch")
precip[,3] <- map[unlist(precip[,3])]

unemp <- unemp %>%
  gather(as.character(states), key = "State", value = "Unemp.Perc")
unemp <- unemp[with(unemp, order(State, Year, Month)),]

pop.density <- pop.density %>%
  gather('2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', key = "Year", value = "Ppl/SqMile")

pop.density$State <- states
pop.density$Year <- as.numeric(pop.density$Year)
pop.density <- pop.density[with(pop.density, order(State, Year)),]

# Cannot find average temperature or precipitation for September, so I'm making them
d <- na.omit(subset(temp, temp$Month == 9))
avg.temp <- aggregate(d[,4], list(d$State), mean)[,2]
var.temp <- diag(aggregate(d[,4], list(d$State), var)[,2])
create.Sept.temp <- round(mvrnorm(1,avg.temp, var.temp),1)
temp$Deg.Far[temp$Year == 2018 & temp$Month == 9] <- create.Sept.temp
temp <- na.omit(temp)

d <- na.omit(subset(precip, precip$Month == 9))
avg.precip <- aggregate(d[,4], list(d$State), mean)[,2]
var.precip <- diag(aggregate(d[,4], list(d$State), var)[,2])
create.Sept.precip <- round(mvrnorm(1,avg.precip, var.precip),1)
precip$Inch[precip$Year == 2018 & precip$Month == 9] <- create.Sept.precip
precip <- na.omit(precip)

# Removing months Jan-Sept for 2010 since we don't have flu data on those monts
temp <- temp[!(temp$Year == 2010 & temp$Month < 10),]
temp <- temp[with(temp, order(State,Year,Month)),]
precip <- precip[!(precip$Year == 2010 & precip$Month < 10),]
precip <- precip[with(precip, order(State,Year,Month)),]

# Repeat the population density values the required number of times so it has the same
# number of rows as temp and precip
pop.density <-data.frame(rep.int(unique(pop.density$State), times=rep(96,47)), rep.int(pop.density$Year, times=rep(c(3, rep(12, 7), 9), 47)), rep.int(pop.density$`Ppl/SqMile`, times=rep(c(3, rep(12, 7), 9), 47)))
colnames(pop.density) <- c("State", "Year", 'Ppl.Sq.Mile')

flu <- read.csv("ILINet_monthly.csv", header = T)
flu[,7] <- flu$ILITotal/flu$TotalPatients
colnames(flu)[7] <- "Proportion"
flu$Proportion[is.na(flu$Proportion)] <- 0


# remove
flu$State = droplevels(flu$State)


# Changing the states into abbrevations
flu$State <- rep(states, each = 12*8)


Y <- flu$Proportion

# Making October (10th month) the baseline
indicators <- matrix(0, nrow = length(flu$Month), ncol = length(unique(flu$Month)))
indicators[cbind(seq_along(flu$Month), flu$Month)] <- 1
indicators <- indicators[,-10]
X <- cbind(1, indicators, temp$Deg.Far, precip$Inch, pop.density$Ppl.Sq.Mile, unemp$Unemp.Perc)

G <- matrix(0, nrow = length(flu$State), ncol = nlevels(flu$State))
G[cbind(seq_along(flu$State), flu$State)] <- 1
