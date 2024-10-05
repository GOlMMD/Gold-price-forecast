#mohammadhosein_Golmohammadi
#replicate time series data of GoldPrice
rm(list=ls(all=TRUE))
ls()
library(readxl)
filepath <- 'D:/Videos/mali/sanji mali/goldprice forcast'
setwd(filepath)
gold <- read_excel("gold.xlsx", col_types = c("text","numeric"))
View(gold)
# load forecast package 
library(lubridate)
library(tseries)
library(timechange)
library(tidyverse)
source("urtests.r")
library(lmtest)
library(stats)

# Export figures? 0 for NO, 1 for YES
exp.plt <- 0
summary(gold)
str(gold)
lGoldPrice<-log(gold$GoldPrice)
lGoldPrice<-as.ts(lGoldPrice)
goldprice<-as.ts(gold$GoldPrice)
gold$Month <- seq(as.Date("2003-11-03"),as.Date("2014-01-14"),length=123,frequency=4)
plot(gold$GoldPrice,col="blue",ylab="GoldPrice",xlab="year")
plot(lGoldPrice,col="blue",ylab="lGoldPrice",xlab="year")
plot(gold)
plot(goldprice)
# plot ACF & PACF
acf(gold$GoldPrice,lag.max=40,main="")
pacf(gold$GoldPrice,lag.max=20,main="")
#AR(1) as initial guess
#checking unit root for stationary
adf.test(lGoldPrice)
adf.test(goldprice)
#Dickey-Fuller = -2.1225, Lag order = 4, p-value = 0.5257
#our data isnt stationary we should get diffrence
#After first diff of our data
laglGoldPrice<-lag(lGoldPrice,-1)
dlGoldPrice<-diff(lGoldPrice)
adf.test(dlGoldPrice)
#Dickey-Fuller = -5.2561, Lag order = 4, p-value = 0.01
#we find out that our data become STATIONARY so integrated part is one
#deterred our real data
trend=lm(gold$GoldPrice~c(1:length(gold$GoldPrice)))
deterend=residuals(trend)
deterend <- ts(deterend,start=c(2003-11-03),frequency=4)
plot.ts(deterend)
test1h <- ur.test(deterend,trend="ct",method="adf.gls",kmax=20)
print.ur.test(test1h)
#that's not even trend stationary
#adf.test of deterend
adf.test(deterend)
acf(deterend,lag.max=40,main="")
pacf(deterend,lag.max=20,main="")
par(mfrow=c(1,2))

# function that computes AIC
aic <- function(model){
  n <- length(model$residuals)
  k <- length(model$coef)
  s2 <- model$sigma2
  return(log(s2)+(2*k/n))
}
# function that computes BIC
bic <- function(model){
  n <- length(model$residuals)
  k <- length(model$coef)
  s2 <- model$sigma2
  rss<-sum(resid(model)^2)
  return(log(rss/n)+(log(n)*k/n))
}
bic2 <- function(model){
  n <- length(model$residuals)
  k <- length(model$coef)
  s2 <- model$sigma2
  rss<-sum(resid(model)^2)
  return(log(s2)+(log(n)*k/n))
}
# max lag orders allowed
p.max <- 3
q.max <- 3

# Storage matrices
aic.results <- matrix(rep(0,(p.max+1)*(q.max+1)),nrow=(q.max+1),ncol=(p.max+1))
bic.results <- matrix(rep(0,(p.max+1)*(q.max+1)),nrow=(q.max+1),ncol=(p.max+1))
bic2.results <- matrix(rep(0,(p.max+1)*(q.max+1)),nrow=(q.max+1),ncol=(p.max+1))
for(i in 0:q.max){
  for(j in 0:p.max){
    model <- arima(gold$GoldPrice,order=c(j,0,i),method="ML")
    aic.results[i+1,j+1] <- aic(model)
    bic.results[i+1,j+1] <- bic(model)
  }
}
bic.results
min(bic.results)
aic.results
min(aic.results)

for(i in 0:q.max){
  for(j in 0:p.max){
    model <- arima(gold$GoldPrice,order=c(j,1,i),method="ML")
    aic.results[i+1,j+1] <- aic(model)
    bic.results[i+1,j+1] <- bic(model)
  }
}
bic.results
min(bic.results)

aic.results
min(aic.results)

# estimate ARima(1,1,1)
# (use functions in the forecast package) 
library(forecast)
model11 <- arima(gold$GoldPrice[1:123],order=c(1,1,1),method="ML")
model11
plot(model11$residuals)
acf(model11$residuals)
pacf(model11$residuals)
fore11 <- forecast(model11,h=6)
layout(1)
plot(fore11)
lines(gold$GoldPrice)
accuracy(fore11)
summary(fore11)
#??????in ja adade dastanesh chie bepors???
polyroot(c(-model11$coef[1],-model11$coef[2],2))
# test for autocorrelation in the residuals  0.5643012+0.2073893i 0.5643012-0.2073893i thats prove stationary
#??????inam kolan dastanesh chie??????
Box.test(model11$resid,lag=18,type="Ljung-Box",fitdf=1)
n.obs <- length(gold$GoldPrice)
n.end <- 117 
n.obs - n.end
# plot sample ACF and PACF 
# (use functions in the forecast package) 
par(mfrow=c(1,2))

Acf(lGoldPrice[1:n.end],lag.max=40)
Pacf(lGoldPrice[1:n.end],lag.max=24)
Acf(gold$GoldPrice[1:n.end],lag.max=40)
Pacf(gold$GoldPrice[1:n.end],lag.max=24)
# estimate AR(2)
# (use functions in the forecast package) 
#thats what we want not  the previous one
model11 <- Arima(gold$GoldPrice[1:n.end],order=c(1,1,1),method="ML")
model11
model12 <- Arima(gold$GoldPrice[1:n.end],order=c(1,1,2),method="ML")
model12
model13 <- Arima(gold$GoldPrice[1:n.end],order=c(1,1,3),method="ML")
model13
model101 <- Arima(gold$GoldPrice[1:n.end],order=c(1,0,1),method="ML")
model101
model102 <- Arima(gold$GoldPrice[1:n.end],order=c(1,0,2),method="ML")
model102
model103 <- Arima(gold$GoldPrice[1:n.end],order=c(1,0,3),method="ML")
model103
# 1-step ahead recursive predictions kolan rideman?????
# set matrix for storage
prede <- matrix(rep(0,6),6,1)
for(i in 1:6){
  x <- goldprice[1:123+i-1] 
  
  model11.tmp <- arima(x,order=c(1,1,1),method="ML")
  prede[i,1] <- forecast(model11.tmp,h=1)$mean[1]
  goldprice[123+i] <- prede[i,1]
}
pred <- matrix(rep(0,36),6,6)
# start loop
for(i in 1:6){
  x <- gold$GoldPrice[1:n.end+i-1] 
  
  model11.tmp <- arima(x,order=c(1,1,1),method="ML")
  pred[i,1] <- forecast(model11.tmp,h=1)$mean[1]
  
  model12.tmp <- arima(x,order=c(1,1,2),method="ML")
  pred[i,2] <- forecast(model12.tmp,h=1)$mean[1]
  
  model13.tmp <- arima(x,order=c(1,1,3),method="ML")
  pred[i,3] <- forecast(model13.tmp,h=1)$mean[1]
  
  model101.tmp <- arima(x,order=c(1,0,1),method="ML")
  pred[i,4] <- forecast(model101.tmp,h=1)$mean[1]
  
  model102.tmp <- arima(x,order=c(1,0,2),method="ML")
  pred[i,5] <- forecast(model102.tmp,h=1)$mean[1]
  
  model103.tmp <- arima(x,order=c(1,0,3),method="ML")
  pred[i,6] <- forecast(model103.tmp,h=1)$mean[1]
}
# set predictions as time series object
pred.ts <- ts(pred,start=c(117),end=c(123))

gold$pred1 <- c(model11$fitted,pred[,1])
gold$pred2 <- c(model12$fitted,pred[,2])
gold$pred3 <- c(model13$fitted,pred[,3])
gold$pred4 <- c(model101$fitted,pred[,4])
gold$pred5 <- c(model102$fitted,pred[,5])
gold$pred6 <- c(model103$fitted,pred[,6])
longer_gold <- pivot_longer(gold, !Month, names_to="model_or_data", values_to = "gold_price")
ggplot(longer_gold[longer_gold$Month>ymd(20040101),]) + # what and when
  geom_line(aes(x=Month, y=gold_price, color=model_or_data)) + # type of graph 
  geom_ribbon(aes(x=Month, y=gold_price, ymin=gold_price-100, ymax=gold_price+100, fill=model_or_data), alpha=0.2) +
  labs(labelsitle = "data vs. model",x = "Date",y = "Price (2003-2014)") +
  scale_x_date(date_breaks = "12 months") +  # x axis ticks
  geom_vline(xintercept = ymd(20130714)) +  # event of interest
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.1, 0.85))
# plot  1-step ahead forecasts
#without ggpolt2
#layout(1)
#plot(gold$GoldPrice)
#lines(pred.ts[,1],col="blue",)
#lines(pred.ts[,2],col="magenta")
#lines(pred.ts[,3],col="red")
#legend(x="topleft",c("ARIma(1,1,3) forecast","ARIma(1,1,1) forecast","ARIMA(1,1,2) forecast"),col=c("red","blue","magenta"),lty=1,lwd=1)

accuracy(pred[,1],gold$GoldPrice[(n.end+1):n.obs])
accuracy(pred[,2],gold$GoldPrice[(n.end+1):n.obs])
accuracy(pred[,3],gold$GoldPrice[(n.end+1):n.obs])
accuracy(pred[,4],gold$GoldPrice[(n.end+1):n.obs])
accuracy(pred[,5],gold$GoldPrice[(n.end+1):n.obs])
accuracy(pred[,6],gold$GoldPrice[(n.end+1):n.obs])
accuracy(model11.tmp)
summary(model11.tmp)
summary(model12.tmp)
summary(model13.tmp)

#pred.ts <- ts(pred,start=c(117),end=c(123))

# plot  1-step ahead forecasts
layout(1)
plot(gold$GoldPrice)
lines(pred.ts[,4],col="blue",)
lines(pred.ts[,5],col="magenta")
lines(pred.ts[,6],col="red")
legend(x="topleft",c("ARIma(1,1,3) forecast","ARIma(1,1,1) forecast","ARIMA(1,1,2) forecast"),col=c("red","blue","magenta"),lty=1,lwd=1)
accuracy(pred[,4],gold$GoldPrice[(n.end+1):n.obs])
accuracy(pred[,5],gold$GoldPrice[(n.end+1):n.obs])
accuracy(pred[,5],gold$GoldPrice[(n.end+1):n.obs])
accuracy(model101.tmp)
summary(model101.tmp)
summary(model102.tmp)
summary(model103.tmp)
modele11 <- Arima(gold$GoldPrice[1:123],order=c(1,1,1),method="ML")
modele11
accuracy(modele11)
modele12 <- Arima(gold$GoldPrice[1:123],order=c(1,1,2),method="ML")
modele12
accuracy(modele12)
modele13 <- Arima(gold$GoldPrice[1:123],order=c(1,1,3),method="ML")
modele13
accuracy(modele13)
modele101 <- Arima(gold$GoldPrice[1:123],order=c(1,0,1),method="ML")
modele101
accuracy(modele101)
modele102 <- Arima(gold$GoldPrice[1:123],order=c(1,0,2),method="ML")
modele102
accuracy(modele102)
modele103 <- Arima(gold$GoldPrice[1:123],order=c(1,0,3),method="ML")
modele103
accuracy(modele103)
Box.test(modele101$resid,lag=18,type="Ljung-Box",fitdf=1)
Box.test(modele102$resid,lag=18,type="Ljung-Box",fitdf=1)
Box.test(modele103$resid,lag=18,type="Ljung-Box",fitdf=1)
Box.test(modele11$resid,lag=18,type="Ljung-Box",fitdf=1)
Box.test(modele12$resid,lag=18,type="Ljung-Box",fitdf=1)
Box.test(modele13$resid,lag=18,type="Ljung-Box",fitdf=1)
sum(goldprice-modele101$fitted)^2/sum((goldprice-mean(goldprice))^2)
sum(goldprice-modele102$fitted)^2/sum((goldprice-mean(goldprice))^2)
sum(goldprice-modele103$fitted)^2/sum((goldprice-mean(goldprice))^2)
sum(goldprice-modele11$fitted)^2/sum((goldprice-mean(goldprice))^2)
sum(goldprice-modele12$fitted)^2/sum((goldprice-mean(goldprice))^2)
sum(goldprice-modele13$fitted)^2/sum((goldprice-mean(goldprice))^2)

