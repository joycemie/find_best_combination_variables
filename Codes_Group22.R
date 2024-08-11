#set up
data = read.csv("C:/Users/陈婕/Desktop/Data_Group22.csv", header = TRUE)
ts<-(data$CUUR0000SEHB)[27:314]
data<-ts(log(ts),start=2000+1/12,frequency=12)
#change to the log
plot(data)

#training and test
train<-window(data, start=c(2000,2), end=c(2022,1))
test<-window(data, start =c(2022,2), end=c(2024,1))
#future time
future.tim<-seq(2024+1/12,2026,by=1/12)
poly.future=poly(future.tim,14,raw=FALSE)
season.future=c(seq(2,12),1,2,seq(3,12),1)

######Regression:

#Unregularized
poly.tim<-poly(time(data),15,raw=FALSE) #orthogonal polynomials
season.train=as.factor(cycle(train)) #seasonal component
apse=c()
#try different degrees
for (i in 1:15) {
  fit<-lm(train~poly.tim[1:264,1:i]+season.train)
  tim.test=poly.tim[265:288,1:i]
  season.test=as.factor(cycle(test))
  model_matrix=model.matrix(test~tim.test+season.test)
  pred=coef(fit)%*%t(model_matrix)
  apse=c(apse,mean((as.vector(test)-pred)^2))
}
#which degree minimize apse
which(apse==min(apse))
#how it fits original plot
tim=time(data)
m.data=as.factor(cycle(data))
model_u<-lm(data~poly.tim[,1:1]+m.data)
plot(data,lwd=1.5, main="Original vs Fitted Data Using Unregularized
     Regression", ylab="log(CPI)", xlim=c(2000,2024.3))
points(time(data),fitted(model_u),type='l',col='blue',lwd=1.5)
legend("bottomright", legend = c("Original","Fitted"), 
       col = c("black","blue"), lty = c(1,1))
#Prediction Using Unregularized Regression
model_u<-lm(data~poly(tim,1)+m.data)
d.frame=data.frame(tim=future.tim,m.data=as.factor(season.future))
pred_u=predict(model_u,d.frame,interval='prediction',level=0.95)
plot(exp(data),lwd=1.5, main="Prediction Using Regression",
     ylab="CPI", xlim=c(2000,2026.3),ylim=c(100,220))
lines(future.tim,exp(pred_u[,'fit']),type='l',col='red',lwd=1.5)
lines(future.tim,exp(pred_u[,'lwr']),type='l',col='blue',lwd=1.5)
lines(future.tim,exp(pred_u[,'upr']),type='l',col='blue',lwd=1.5)
legend("bottomright", legend = c("Original","Predict","Interval"), 
       col = c("black","red","blue"), lty=c(1,1,1))

#Residuals
u_res=model_u$residuals
par(mfrow = c(1,2))
plot(u_res,type="p",main="Residuals for Unregularized Regression")
acf(u_res,main="ACF for Residuals")
#try degree 3 since apse does not have a big difference
model_u<-lm(data~poly.tim[,1:3]+m.data)
par(mfrow = c(1,1))
plot(data,lwd=1.5, main="Original vs Fitted Data Using Unregularized
     Regression", ylab="log(CPI)", xlim=c(2000,2025))
points(time(data),fitted(model_u),type='l',col='blue',lwd=1.5)
u_res=model_u$residuals
par(mfrow = c(1,2))
plot(u_res,type="p",main="Residuals for Unregularized Regression")
acf(u_res,main="ACF for Residuals")
#try differencing
diff_1=diff(u_res)
plot(diff_1,main="Residuals after Differencing")
acf(diff_1,"ACF of Residuals")

#Regularized
library(glmnet)
#alpha=0 ridge
lambda_p_rid=c()
pred_mse_rid=c()
for (i in 1:15){
  set.seed(1)
  train.index=as.matrix(cbind(poly.tim[1:264,1:i],season.train))
  fit_rid=cv.glmnet(train.index,as.matrix(train),alpha=0,
                    nfolds=10,type='mse')
  lambda_p_rid<-c(lambda_p_rid,fit_rid$lambda.1se)
  j=which(fit_rid$lambda.1se==fit_rid$lambda)
  mse=fit_rid$cvm[j]
  pred_mse_rid<-c(pred_mse_rid,mse)
}
#which model is the best
which(pred_mse_rid==min(pred_mse_rid))
#fit the model on data
tim.index=as.matrix(cbind(poly.tim[,1:15],m.data))
fit.ridge=glmnet(tim.index, data, alpha=0, 
                 lambda=lambda_p_rid[15], standardize=TRUE)
rid_pred=predict(fit_rid,newx=tim.index)
par(mfrow = c(1,1))
plot(data, type="l", col="black", lwd=1.5, 
     ylab="log(CPI)", xlab = "Time",
     main="Original vs. Fitted Data, alpha=0")
points(time(data),rid_pred,type='l',col='red',lwd = 1.5)
legend("bottomright", legend = c("Original","Fitted"), 
       col = c("black","red"), lty = c(1, 1))

#alpha=1 Lasso
lambda_p_las=c()
pred_mse_las=c()
for (i in 1:15){
  set.seed(1)
  train.index=as.matrix(cbind(poly.tim[1:264, 1:i], season.train))
  fit_las=cv.glmnet(train.index,as.matrix(train),alpha=1,
                    nfolds=10,type='mse')
  lambda_p_las<-c(lambda_p_las,fit_las$lambda.1se)
  j=which(fit_las$lambda.1se==fit_las$lambda)
  mse=fit_las$cvm[j]
  pred_mse_las<-c(pred_mse_las,mse)
}
#which model is the best
which(pred_mse_las==min(pred_mse_las))
#fit the model on data
tim.index=as.matrix(cbind(poly.tim[,1:14],m.data))
fit_las=glmnet(tim.index, data, alpha=1,
               lambda=lambda_p_las[14], standardize=TRUE)
las_pred=predict(fit_las,newx=as.matrix(tim.index))
par(mfrow = c(1,1))
plot(data, type="l", col="black", lwd=1.5, 
     ylab="log(CPI)", xlab = "Time",
     main="Original vs. Fitted Data, alpha=1")
points(time(data),las_pred,type='l',col='red',lwd = 1.5)
legend("bottomright", legend = c("Original","Fitted"), 
       col = c("black","red"), lty = c(1, 1))


#alpha=0.5 Elastic
lambda_p_ela=c()
pred_mse_ela=c()
for (i in 1:15){
  set.seed(1)
  train.index=as.matrix(cbind(poly.tim[1:264, 1:i], season.train))
  fit_ela=cv.glmnet(train.index,as.matrix(train),alpha=0.5,
                    nfolds=10,type='mse')
  lambda_p_ela<-c(lambda_p_ela,fit_ela$lambda.1se)
  j=which(fit_ela$lambda.1se==fit_ela$lambda)
  mse=fit_ela$cvm[j]
  pred_mse_ela<-c(pred_mse_ela,mse)
}
#which model is the best
which(pred_mse_ela==min(pred_mse_ela))
#fit the model on data
tim.index=as.matrix(cbind(poly.tim[,1:12],m.data))
fit_ela=glmnet(tim.index,data,alpha=0.5,
               lambda=lambda_p_ela[12], standardize=TRUE)
ela_pred=predict(fit_ela,newx=as.matrix(tim.index))
par(mfrow = c(1,1))
plot(data, type="l", col="black", lwd=1.5, 
     ylab="CPI", xlab = "Time",
     main="Original vs. Fitted Data, alpha=0.5")
points(time(data),ela_pred,type='l',col='red',lwd = 1.5)
legend("bottomright", legend = c("Original","Fitted"), 
       col = c("black","red"), lty = c(1, 1))

# APSE
test.index=cbind(poly.tim[265:288,1:15],season.test)
pred_ridge=predict(fit_rid,newx=test.index,type="response")
mean((test-pred_ridge)^2) #regularized(ridge)

test.index=cbind(poly.tim[265:288,1:14],season.test)
pred_lasso=predict(fit_las,newx=test.index,type="response")
mean((test-pred_lasso)^2) # regularized(lasso)

test.index=cbind(poly.tim[265:288,1:12],season.test)
pred_elastic=predict(fit_ela,newx=test.index,type="response")
mean((test-pred_elastic)^2) # regularized(elastic)

# Prediction Using Lasso
future.index=cbind(poly.future[,1:14],as.factor(season.future))
pred_las=predict(fit_las,future.index)
plot(exp(data), type="l", col="black", lwd=1.5, ylim=c(75,270),
     ylab="CPI", xlab = "Time", xlim=c(2000,2026.3),
     main="Prediction Using Lasso")
lines(future.tim,exp(pred_las),lwd=1.5,col='red')
legend("topleft", legend = c("Original","Predicted"), 
        col = c("black","red"), lty = c(1, 1))

# Prediction Using Elastic net
future.index=cbind(poly.future[,1:12],as.factor(season.future))
pred_ela=predict(fit_ela,future.index)
plot(exp(data), type="l", col="black", lwd=1.5, ylim=c(75,335),
     ylab="CPI", xlab = "Time", xlim=c(2000,2026.3),
     main="Prediction Using Elastic Net")
lines(future.tim,exp(pred_ela),lwd=1.5,col='red')
legend("topleft", legend = c("Original","Predicted"), 
       col = c("black","red"), lty = c(1, 1))


######residuals
#Ridge
rid_res=data-rid_pred
plot(rid_res,type="l",main="Residuals for Ridge Regression")
acf(rid_res,lag.max=70,main="ACF for Residuals (Ridge)")
#Lasso
las_res=data-las_pred
plot(las_res,type="l",main="Residuals for Lasso Regression")
acf(las_res,main="ACF for Residuals (Lasso)")
#Elastic
ela_res=data-ela_pred
plot(ela_res,type='l',main="Residuals for Elastic Regression")
acf(ela_res,main="ACF for Residuals (Elastic)")

#try differencing
par(mfrow = c(1,1))
diff_one=diff(ela_res,lag=12)
plot(diff_one)
acf(diff_one,lag.max = 24)
diff_two=diff(diff_one,differences=1)
plot(diff_two)
acf(diff_two,lag.max=24)


#Smoothing:

# Regular Differencing
diff_data <- diff(data, differences = 1)
# Seasonal and trend Differencing
seasonal_diff_data <- diff(diff_data, lag = 12)
# Plot the differenced data
par(mfrow = c(1, 2)) # Set up the plotting area to display 2 plots vertically
plot(diff_data, main = "Regular Differenced Data")
plot(seasonal_diff_data, main = "Seasonally Differenced Data")

par(mfrow = c(1, 1))
# ACF plot to check for stationarity
acf(seasonal_diff_data, lag.max=12*24, 
    main = "ACF of Seasonality and Trend Differenced Data")

par(mfrow = c(1, 2))
# Perform additive decomposition
additive_decomp <- decompose(data ,  type = "additive")
# Plot the additive decomposition components
plot(additive_decomp)

# Perform multiplicative decomposition
multiplicative_decomp <- decompose(data, type = "multiplicative")
# Plot the multiplicative decomposition components
plot(multiplicative_decomp)

#Since it has the trend and seasonality, so it could not be use to estimate or prediction.

#HW method
hw.additive = HoltWinters(train, seasonal = "additive")
hw.predict = predict(hw.additive, n.ahead = 24)
mse1 = mean((test - hw.predict)^2)
mse1

hw.mult = HoltWinters(train, seasonal = "multiplicative")
hw.predict = predict(hw.mult, n.ahead = 24)
mse2 = mean((test - hw.predict)^2)
mse2

par(mfrow = c(1,1))
hw.final = HoltWinters(data, seasonal = "multiplicative")
plot(hw.final, ylim = c(4.5, 5.4), 
     predict(hw.final, n.ahead = 12,  prediction.interval = TRUE))


#BJ:
data = read.csv("Data_Group22.csv", header = TRUE)
data = data[-c(1:25),] # delete data before 2000
data = ts(data[-1], start = 1999 + 12/12, frequency = 12)

training <- window(data,2000,2022)
testing <- window(data, 2022+1/12)

plot(log(training), main = "raw log(training data)")
acf(log(training)) # not stationary, trend and seasonality
pacf(log(training)) # exponential decay 


# Seasonal Differencing:
delta_12 <- diff(log(training), lag = 12)
plot(delta_12, main = "Seasonal Differencing")
acf(delta_12, main = "Seasonal Differencing") 
pacf(delta_12,main = "Seasonal Differencing")

# Regular differencing (remove trend)
regu_delta_12 <- diff(delta_12, differences = 1)
plot(regu_delta_12, main = "Regular Differencing")
acf(regu_delta_12, main = "Regular Differencing") 
pacf(regu_delta_12,main = "Regular Differencing")

library(astsa)
# After regular and seasonal differencing, we established d = 1, D = 1, s = 12
# Then, 3 Candidate models: 
# SARIMA(0,1,0) × (1,1,1) & SARIMA(0,1,0) × (0,1,1) & SARIMA(0,1,0) × (2,1,0)

fit1 <- sarima(log(training), p=0,d=1,q=0,P=1,D=1,Q=1,S=12)

fit2 <- sarima(log(training), p=0,d=1,q=0,P=0,D=1,Q=1,S=12)

fit3 <- sarima(log(training), p=0,d=1,q=0,P=2,D=1,Q=0,S=12)

a <- rbind(fit1$ICs, fit2$ICs, fit3$ICs)
data.frame(row.names = c("SARIMA(0,1,0) × (1,1,1)_12",
                         "SARIMA(0,1,0) × (0,1,1)_12", 
                         "SARIMA(0,1,0) × (2,1,1)_12"), 
           AIC = a[,1], AICc = a[,2], BIC = a[,3])

# Forecasting
par(mfrow = c(3,1))
fore1 <- sarima.for(log(training), n.ahead=24, 
                    p=0,d=1,q=0,P=1,D=1,Q=1,S=12)
title("SARIMA(0,1,0)x(1,1,1)_12")

fore2 <- sarima.for(log(training), n.ahead=24, 
                    p=0,d=1,q=0,P=0,D=1,Q=1,S=12)
title("SARIMA(0,1,0)x(0,1,1)_12")

fore3 <- sarima.for(log(training), n.ahead=24, 
                    p=0,d=1,q=0,P=2,D=1,Q=0,S=12)
title("SARIMA(0,1,0)x(2,1,0)_12")

# forecasting based on SARIMA(0,1,0)x(1,1,1)_12
lower1 <- fore1$pred-1.96*fore1$se
upper1 <- fore1$pred+1.96*fore1$se
fit1 <- fore1$pred
plot(data,xlim=c(2000,2025),ylim=c(100,240), main='SARIMA(0,1,0)x(1,1,1)_12')
lines(exp(fit1),col='red',type='b',pch='*') 
lines(exp(lower1),col='blue',lty=2)
lines(exp(upper1),col='blue',lty=2)

lower2 <- fore2$pred-1.96*fore2$se
upper2 <- fore2$pred+1.96*fore2$se
fit2 <- fore2$pred
plot(data,xlim=c(2000,2025),ylim=c(100,240), main='SARIMA(0,1,0)x(0,1,1)_12')
lines(exp(fit2),col='red',type='b',pch='*') 
lines(exp(lower2),col='blue',lty=2)
lines(exp(upper2),col='blue',lty=2)

lower3 <- fore3$pred-1.96*fore3$se
upper3 <- fore3$pred+1.96*fore3$se
fit3 <- fore3$pred
plot(data,xlim=c(2000,2025),ylim=c(100,240), main='SARIMA(0,1,0)x(2,1,0)_12')
lines(exp(fit3),col='red',type='b',pch='*') 
lines(exp(lower3),col='blue',lty=2)
lines(exp(upper3),col='blue',lty=2)

f1 <- mean((fore1$pred-log(testing))^2)
f2 <- mean((fore2$pred-log(testing))^2)
f3 <- mean((fore3$pred-log(testing))^2)

data.frame(row.names = c("SARIMA(0,1,0) × (1,1,1)_12",
                         "SARIMA(0,1,0) × (0,1,1)_12", 
                         "SARIMA(0,1,0) × (2,1,1)_12"), APSE = c(f1,f2,f3))

# our final model is SARIMA(0,1,0) × (1,1,1)_12

# 2 year forecasting

FinalFit <- sarima(log(data),p=0,d=1,q=0,P=1,D=1,Q=1,S=12)
future.forecast <- sarima.for(log(data), n.ahead=24, 
                              p=0,d=1,q=0,P=1,D=1,Q=1,S=12)
lower <- future.forecast$pred-1.96*future.forecast$se
upper <- future.forecast$pred+1.96*future.forecast$se
fit <- future.forecast$pred
ts.plot(data,xlim=c(2000,2026),ylim=c(100, 245),
        main='SARIMA(0,1,0)x(1,1,1)_12')

x = c(time(upper) , rev(time(upper)))
y = c(exp(upper) , rev(exp(lower)))
polygon(x, y, col = "grey" , border =NA)

lines(exp(fit),col='red',type='b',pch=16 , cex=0.5) 
lines(exp(lower),col='black',lty=2)
lines(exp(upper),col='black',lty=2)

#conclusion
data.frame(Name = c("Regression", "Smoothing", "Box-Jenkins"),
           Type = c("Elastic net", "Multiplicate", "SARIMA(0,1,0) × (1,1,1)_12"),
           APSE = c(0.002841686, 0.002847429, 0.001051204))
