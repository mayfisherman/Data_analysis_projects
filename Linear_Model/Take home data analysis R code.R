## load data and extract the response and the design matrix
health.data <- read.table("D:\\Course material\\Statistics and population genetics\\STA 721\\Take home Data analysis\\costs.txt", header = T)
names(health.data)
plot(health.data)
Y <- health.data[,1]
X <- data.matrix(health.data[,2:8])
X.norm <- scale(X)

## compute the eigenvalues of X^T X, where X here is the full design matrix including the intercept
X.extend <- cbind(1,X)
eigen(t(X.extend) %*% X.extend)$values

## Frequentist normal regression and diagnostic plots
cost.lm <- lm(COST ~ RXPM + GS + RI + COPAY + AGE + F + MM, data = health.data)
summary(cost.lm)
par(mfrow=c(2,2))
plot(cost.lm)

## correlation between variables
health.cor <- cor(health.data[-9])
health.cor
library(lattice)
levelplot(health.cor)

## Zellner's g prior with g = n using the full model
library(BMS)
cost.zlm <- zlm(COST ~ RXPM + GS + RI + COPAY + AGE + F + MM, data = health.data, g = "UIP")
summary(cost.zlm)
coef(cost.zlm)
## added-variable plots for g-prior regression
library(car)
avPlots(cost.zlm, main="Added-Variable Plot: Full model with g prior")

## ridge regression 
library(MASS)
cost.ridge = lm.ridge(Y ~ X.norm, lambda = seq(0, 10, 0.0001))
best.lambda = as.numeric(names(which.min(cost.ridge$GCV)))
best.lambda
cost.ridge.refit = lm.ridge(Y ~ X.norm, lambda = best.lambda)
summary(cost.ridge.refit)
coef(cost.ridge.refit)
plot(cost.ridge.refit)

## lasso regression
library(lars)
cost.lasso <- lars(X.norm, Y, type="lasso")
cost.Cp <- summary(cost.lasso)$Cp
best.Cp <- (1:length(cost.Cp))[cost.Cp == min(cost.Cp)]
coef.lasso = coef(cost.lasso)[best.Cp,]
coef.lasso 
plot(cost.lasso)

## horseshoe regression
library(monomvn)
cost.bla <- blasso(X.norm, Y, T=11000, case = "hs", RJ = FALSE, normalize = F)
plot(cost.bla, burnin = 1000)
## give posterior summaries using post-burnin MCMC samples
coef.bla = c(mean(cost.bla$mu[1001:11000]),apply(cost.bla$beta[1001:11000,], 2, mean))
quantile(cost.bla$mu[1001:11000],c(.025,.975))
quantile(cost.bla$beta[1001:11000,1],c(.025,.975))
quantile(cost.bla$beta[1001:11000,2],c(.025,.975))
quantile(cost.bla$beta[1001:11000,3],c(.025,.975))
quantile(cost.bla$beta[1001:11000,4],c(.025,.975))
quantile(cost.bla$beta[1001:11000,5],c(.025,.975))
quantile(cost.bla$beta[1001:11000,6],c(.025,.975))
quantile(cost.bla$beta[1001:11000,7],c(.025,.975))

## Bayesian model averaging with the Zellner-Siow Cauchy prior
library(BAS)
cost.bma <- bas.lm(COST ~ RXPM + GS + RI + COPAY + AGE + F + MM, data = health.data, prior="ZS-null", alpha=3, n.models=2^7, update=50, initprobs="eplogp")
par(mfrow=c(2,2))
plot(cost.bma, ask = F)
par(mfrow=c(1,1))
image(cost.bma)
coef.bma = coef(cost.bma)
coef.bma
par(mfrow=c(2,4))
plot(coef.bma, subset=1:8,ask=F)
summary(cost.bma)  ## show results with the top 5 models




## load data and extract the response and the design matrix
health.data <- read.table("D:\\Course material\\Statistics and population genetics\\STA 721\\Take home Data analysis\\costs.txt", header = T)
names(health.data)
plot(health.data)
Y <- health.data[,1]
X <- data.matrix(health.data[,2:8])
X.norm <- scale(X)

## log transform the MM in original data
attach(health.data)
health.data$logMM <- log(MM) 
detach(health.data)
health.newdata = data.frame(health.data[,c(1:7, 10)])

plot(health.newdata)
Y.new <- health.newdata[,1]
X.new <- data.matrix(health.newdata[,2:8])
X.newnorm <- scale(X.new)

## compute the eigenvalues of X^T X, where X here is the full design matrix including the intercept
X.newextend <- cbind(1,X.new)
eigen(t(X.newextend) %*% X.newextend)$values

## Frequentist normal regression and diagnostic plots
cost.newlm <- lm(COST ~ ., data = health.newdata)
summary(cost.newlm)
par(mfrow=c(2,2))
plot(cost.newlm)
vif(cost.newlm)

## correlation between variables
health.newcor <- cor(X.new)
health.newcor
library(lattice)
levelplot(health.newcor)

## Zellner's g prior with g = n using the full model
library(BMS)
cost.newzlm <- zlm(COST ~ ., data = health.newdata, g = "UIP")
summary(cost.newzlm)
coef(cost.newzlm)
## added-variable plots for g-prior regression
library(car)
avPlots(cost.newzlm, main="Added-Variable Plot: Full model with g prior")

## ridge regression 
library(MASS)
cost.newridge = lm.ridge(Y.new ~ X.newnorm, lambda = seq(0, 10, 0.0001))
best.newlambda = as.numeric(names(which.min(cost.newridge$GCV)))
best.newlambda
cost.ridge.newrefit = lm.ridge(Y.new ~ X.newnorm, lambda = best.newlambda)
summary(cost.ridge.newrefit)
coef(cost.ridge.newrefit)
plot(cost.ridge.newrefit)

## lasso regression
library(lars)
cost.newlasso <- lars(X.newnorm, Y.new, type="lasso")
cost.newCp <- summary(cost.newlasso)$Cp
best.newCp <- (1:length(cost.newCp))[cost.newCp == min(cost.newCp)]
coef.newlasso = coef(cost.newlasso)[best.newCp,]
coef.newlasso 
plot(cost.newlasso)

## horseshoe regression
library(monomvn)
cost.newbla <- blasso(X.newnorm, Y.new, T=11000, case = "hs", RJ = FALSE, normalize = F)
plot(cost.newbla, burnin = 1000)
## give posterior summaries using post-burnin MCMC samples
coef.newbla = c(mean(cost.newbla$mu[1001:11000]),apply(cost.newbla$beta[1001:11000,], 2, mean))
coef.newbla
quantile(cost.newbla$mu[1001:11000],c(.025,.975))
quantile(cost.newbla$beta[1001:11000,1],c(.025,.975))
quantile(cost.newbla$beta[1001:11000,2],c(.025,.975))
quantile(cost.newbla$beta[1001:11000,3],c(.025,.975))
quantile(cost.newbla$beta[1001:11000,4],c(.025,.975))
quantile(cost.newbla$beta[1001:11000,5],c(.025,.975))
quantile(cost.newbla$beta[1001:11000,6],c(.025,.975))
quantile(cost.newbla$beta[1001:11000,7],c(.025,.975))

## Bayesian model averaging with the Zellner-Siow Cauchy prior
library(BAS)
cost.newbma <- bas.lm(COST ~ ., data = health.newdata, prior="ZS-null", alpha=3, n.models=2^7, update=50, initprobs="eplogp")
par(mfrow=c(2,2))
plot(cost.newbma, ask = F)
par(mfrow=c(1,1))
image(cost.newbma)
coef.newbma = coef(cost.newbma)
coef.newbma
par(mfrow=c(2,4))
plot(coef.newbma, subset=2:8,ask=F)
summary(cost.newbma)  ## show results with the top 5 models

