#
########################################################################
#                                                                      #
#                         EUR_USD t-dist                               #
#                                                                      #
########################################################################
#
##
## Data Prep
rm(list=ls())
setwd('D:\\RICE UNIVERSITY\\621\\ASSIGNMENT\\6')
library(FinTS); library(tseries)
library(rugarch); library(fGarch)
library(forecast); library(stochvol)
library(lubridate);library(parallel)
library(doParallel);library(foreach)
EUR_USD_all <- read.csv("EUR_USD.csv", header = TRUE)
EUR_USD_all <- EUR_USD_all[!is.na(EUR_USD_all$VALUE),]
EUR_USD <- EUR_USD_all[1:2764, ]
EUR_USD_pred <- EUR_USD_all[-c(1:2764), ]
## t distribution error
ret.all.t <- logret(EUR_USD_all$VALUE, demean = TRUE)
ret.t <- ret.all.t[1:2763]
res.t <- svsample(ret.t,designmatrix = "ar0", priormu = c(-10, 1),
                  priorphi = c(20, 1.1),
                  priorsigma = .1,
                  priornu = c(2,100))
summary(res.t, showlatent = FALSE)
# volplot(res.t, forecast = 100, dates = EUR_USD_all$DATE[-1])
plot(res.t, showobs = FALSE)
myresid.t <- resid(res.t)
plot(myresid.t, ret.t)
## Prediction
ret.pred.t <- ret.all.t[-c(1:2763)]
PI.t <- data.frame(Date = EUR_USD_pred$DATE,
                   lower = rep(0, times = length(ret.pred.t)),
                   upper = rep(0, times = length(ret.pred.t)))
ret.new.t <- ret.t
for(i in 1:length(ret.pred.t)) {
  predvol.t <- predict(res.t, 1)
  preddraws.t <- arpredict(res.t, predvol.t)
  predquants.t <- apply(preddraws.t, 2, quantile, c(.1,.9))
  PI.t$lower[i] <- predquants.t[1]
  PI.t$upper[i] <- predquants.t[2] ## Prediction interval
  ret.new.t <- c(ret.new.t,ret.pred.t[i])
  res.t <- svsample(ret.new.t,designmatrix = "ar0", 
                    priormu = c(-10, 1),
                    priorphi = c(20, 1.1), ## c(0.5, 0.5)
                    priorsigma = .1,
                    priornu = c(2,100))
}
plot.PI.t <- PI.t
plot.PI.t$real <- ret.pred.t
write.csv(plot.PI.t, file = "plot.PI.t.csv")
cover.t <- (length(which(plot.PI.t$real>=PI.t$lower & 
                           plot.PI.t$real<=PI.t$upper)))/
  nrow(plot.PI.t)
## Visualization
plot(real ~ Date, plot.PI.t, main = "Observed and predicted log return 
     for t distribution error", 
     col = 3)
lines(plot.PI.t$lower, col = 4, lty = 2)
lines(plot.PI.t$upper, col = 4, lty = 2)
##
############## sensitivity check ##########
##
##
## Data Prep
rm(list=ls())
EUR_USD_all <- read.csv("EUR_USD.csv", header = TRUE)
EUR_USD_all <- EUR_USD_all[!is.na(EUR_USD_all$VALUE),]
EUR_USD <- EUR_USD_all[1:2764, ]
EUR_USD_pred <- EUR_USD_all[-c(1:2764), ]
## change priorphi
## t distribution error
ret.all.t2 <- logret(EUR_USD_all$VALUE, demean = TRUE)
ret.t2 <- ret.all.t2[1:2763]
res.t2 <- svsample(ret.t2,designmatrix = "ar0", 
                   priormu = c(-10, 1),
                   priorphi = c(0.5, 0.5),
                   priorsigma = .1,
                   priornu = c(2,100))
summary(res.t2, showlatent = FALSE)
# volplot(res.t, forecast = 100, dates = EUR_USD_all$DATE[-1])
plot(res.t2, showobs = FALSE)
myresid.t2 <- resid(res.t2)
plot(myresid.t2, ret.t2)
## Prediction
ret.pred.t2 <- ret.all.t2[-c(1:2763)]
PI.t2 <- data.frame(Date = EUR_USD_pred$DATE,
                    lower = rep(0, times = length(ret.pred.t2)),
                    upper = rep(0, times = length(ret.pred.t2)))
ret.new.t2 <- ret.t2
for(i in 1:length(ret.pred.t2)) {
  predvol.t2 <- predict(res.t2, 1)
  preddraws.t2 <- arpredict(res.t2, predvol.t2)
  predquants.t2 <- apply(preddraws.t2, 2, quantile, c(.1,.9))
  PI.t2$lower[i] <- predquants.t2[1]
  PI.t2$upper[i] <- predquants.t2[2] ## Prediction interval
  ret.new.t2 <- c(ret.new.t2,ret.pred.t2[i])
  res.t2 <- svsample(ret.new.t2,designmatrix = "ar0", 
                     priormu = c(-10, 1),
                     priorphi = c(0.5, 0.5),
                     priorsigma = .1,
                     priornu = c(2,100))
}
plot.PI.t2 <- PI.t2
plot.PI.t2$real <- ret.pred.t2
write.csv(plot.PI.t2, file = "plot.PI.t2.csv")
cover.t2 <- (length(which(plot.PI.t2$real>=PI.t2$lower & 
                            plot.PI.t2$real<=PI.t2$upper)))/
  nrow(plot.PI.t2)
## Visualization
plot(real ~ Date, plot.PI.t2, main = "Observed and predicted log return 
     for t distribution error2", 
     col = 3)
lines(plot.PI.t2$lower, col = 4, lty = 2)
lines(plot.PI.t2$upper, col = 4, lty = 2)