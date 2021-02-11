library(maptools)
library(spdep)
library(maps)

##############################################################
##############################################################
#
# Neighborhood matrix for states
#
usa.state = map(database="state", fill=TRUE, plot=FALSE)
state.ID <- sapply(strsplit(usa.state$names, ":"), function(x) x[1])
usa.poly = map2SpatialPolygons(usa.state, IDs=state.ID)
usa.nb = poly2nb(usa.poly)
usa.adj.mat = nb2mat(usa.nb, style="B")

##Write the 0-1 adjacency matrix
W.raw = usa.adj.mat
#Remove D.C. and Florida
W=W.raw[c(1:7,10:49),c(1:7,10:49)]
colnames(W) <- rownames(W)
D = diag(colSums(W))


######################################################
######################################################
#
# Reading in our data
#

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




#################################################################################
# Necessary Packages

library(Matrix) 
library(hierarchicalDS)
library(mvtnorm)
library(mnormt)
library(lattice)


#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
# 
# Full block Gibbs sampler:
# 
#


flu.mcmc<-function(Y,X,D,W,rho,b,sig2,tau2,as2,bs2,at2,bt2,N,G,verbose=FALSE){
  
  V<-dim(D)[1]
  n <- length(Y)
  DW<-D-rho*W
  I<-sparseMatrix(1:V,1:V,x=1) 
  Ri <- solve(R)
  
  b.mcmc<-matrix(-99,nrow=N,ncol=V)
  sig2.mcmc<-rep(-99,N)
  tau2.mcmc<-rep(-99,N)
  beta0.mcmc<-matrix(-99,nrow = N,ncol=dim(X)[2])
  
  
  for(g in 1:N){
    
    ##############################################
    # Sample the intercept
    
    XtXR <- solve(t(X)%*%X + sig2*Ri)
    Y.s<-Y-G%*%b
    
    beta0<-mvrnorm(1, XtXR%*%t(X)%*%Y.s, sig2*XtXR)
    
    ##############################################
    # sample spatial random effects 
    
    Yt<-t(G)%*%(Y-X%*%beta0)
    
    Prec<- t(G)%*%G + sig2*DW/(tau2)
    CH<-chol(Prec)
    R1<-solve(CH,rnorm(V,0,sqrt(sig2)),sparse=TRUE)
    R2<-solve(t(CH),Yt,sparse=TRUE)
    R3<-solve(CH,R2,sparse=TRUE)
    b<-as.vector(R1+R3)
    
    
    
    ###############################################
    # Sample sig2
    
    as2s<-as2+n/2
    bs2s<-as.vector(bs2+ sum((Y-X%*%beta0-G%*%b)^2)/2)
    
    samp<-rgamma(1,as2s,bs2s)
    
    sig2<-1/samp
    
    ###############################################
    # Sample tau2
    
    at2s<-at2+V/2
    bt2s<-as.vector(bt2+ t(b)%*%DW%*%b/(2))
    
    samp<-rgamma(1,at2s,bt2s)
    
    tau2<-1/samp
    
    sig2.mcmc[g]<-sig2
    tau2.mcmc[g]<-tau2
    beta0.mcmc[g,]<-beta0
    b.mcmc[g,]<-b
    
    if(verbose==TRUE){
      print(g)
      
      if(g%%100==0){
        par(mfrow=c(2,1))
        plot(sig2.mcmc[1:g])
        plot(tau2.mcmc[1:g])
      }
    } 
  }
  
  return(list("sig2"=sig2.mcmc,"tau2"=tau2.mcmc,"beta0"=beta0.mcmc,"b"=b.mcmc))
}




#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
# 
# Initialize parameters
# 
#

# Number of Iterations
N = 10000

# Beta
R = diag(rep(100,ncol(X)))

# Spatial
rho = 0.9
b = rep(0,47)

# Sigma2
sig2 = 2
as2 = 0
bs2 = 0


# Tau2
tau2 = 2
at2 = 0
bt2 = 0



#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
# 
# Run the model
# 
#
flu.gibbs <- flu.mcmc(Y,X,D,W,rho,b,sig2,tau2,as2,bs2,at2,bt2,N,G,verbose=FALSE)

# Trace plot of beta's
par(mfrow = c(4,4))
for(i in 1:ncol(flu.gibbs$beta0)){
  plot(flu.gibbs$beta0[,i], ylab = paste("Beta",as.character(i-1)), main = paste("Trace Plot of Beta",i-1))
}

# Means of beta's
colMeans(flu.gibbs$beta0)

#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
# 
# Plotting chosen states' predicted and true values, one from each region of US
# 
#
par(mfrow = c(2,3))
ch.states = c('SC', 'PA', 'TX', 'SD', 'UT', 'ID')
for(i in 1:length(ch.states)){
  hold.x <- data.frame(flu$State, X)
  hold.x <- subset(hold.x, hold.x[,1] == as.character(ch.states[i]))
  hold.y <- subset(flu, flu$State == ch.states[i])
  plot(hold.y$Proportion, xlab = "Months since Oct 2010", ylab = "Proportion of Ppl with the Flu", main = paste("True and Predicted values for", as.character(ch.states[i])), type = 'l')
  points(1:96, as.matrix(hold.x[,-1])%*%colMeans(flu.gibbs$beta0) + colMeans(flu.gibbs$b)[which(states == as.character(ch.states[i]))], col = 'red')
}


HPDinterval(as.mcmc(flu.gibbs$beta0), cred = 0.95) # var 6,14,15 not significant
HPDinterval(as.mcmc(flu.gibbs$sig2), cred = 0.95)
HPDinterval(as.mcmc(flu.gibbs$tau2), cred = 0.95)


##############################################################################################
##############################################################################################
##############################################################################################
#
# Training/Testing
#

X.train <- data.frame(flu$Year, X)
X.train <- as.matrix(subset(X.train, X.train$flu.Year != 2018 & !( flu$Year == 2017 & flu$Month %in% c(10,11,12)))[,-1])
Y.train <- flu[flu$Year != 2018 & !( flu$Year == 2017 & flu$Month %in% c(10,11,12)),]$Proportion
G <- matrix(0, nrow = length(Y.train), ncol = nlevels(flu$State))
G[cbind(seq_along(Y.train), flu$State[1:length(Y.train)])] <- 1

flu.gibbs.train <- flu.mcmc(Y.train,X.train,D,W,rho,b,sig2,tau2,as2,bs2,at2,bt2,N,G,verbose=FALSE)

X.test <- data.frame(flu$Year, X)
X.test <- subset(X.test, X.test$flu.Year == 2018 | (flu$Year == 2017 & flu$Month %in% c(10,11,12)))[,-1]

par(mfrow = c(2,3))
ch.states = c('SC', 'PA', 'TX', 'SD', 'UT', 'ID')
for(i in 1:length(ch.states)){
  hold.x <- data.frame(rep.int(states,times = rep(12,47)), X.test)
  hold.x <- subset(hold.x, hold.x[,1] == as.character(ch.states[i]))
  hold.y <- subset(flu, flu$State == ch.states[i] & (flu$Year == 2018 | (flu$Year == 2017 & flu$Month %in% c(10,11,12))))
  plot(hold.y$Proportion, xlab = "Months since Oct 2018", ylab = "Proportion of Ppl with the Flu", main = paste("True and Predicted values for", as.character(ch.states[i])), type = 'l')
  points(1:12, as.matrix(hold.x[,-1])%*%colMeans(flu.gibbs.train$beta0) + colMeans(flu.gibbs.train$b)[which(states == as.character(ch.states[i]))], col = 'red')
}
