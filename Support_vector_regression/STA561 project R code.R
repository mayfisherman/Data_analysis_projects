library(MASS)
library(lars)
library(glmnet)
library(e1071)

setwd("D:/Course material/Statistics and population genetics/STA 561/project/data")
house.data <- read.table('housing_data.txt', header = F)
colnames(house.data) <- c('CRIM', 'ZN', 'INDUS', 'CHAS', 'NOX', 'RM', 'AGE', 'DIS', 'RAD', 'TAX', 'PTRATIO', 'B', 'LSTAT', 'MEDV')
house.lm1 <- lm(MEDV ~ ., data = house.data)

MSE.svm <- rep(NULL, 10)
R2.svm <- rep(NULL, 10)
MSE.svml <- rep(NULL, 10)
R2.svml <- rep(NULL, 10)
MSE.mlr <- rep(NULL, 10)
R2.mlr <- rep(NULL, 10)
MSE.ridge <- rep(NULL, 10)
R2.ridge <- rep(NULL, 10)
MSE.lasso <- rep(NULL, 10)
R2.lasso <- rep(NULL, 10)
MSE.lassomin <- rep(NULL, 10)
R2.lassomin <- rep(NULL, 10)


for (i in 1:10){
  index <- 1:nrow(house.data)
  testindex <- sample(index, trunc(length(index)/3))
  testset <- house.data[testindex,]
  trainset <- house.data[-testindex,]
  
  # SVR using a RBF kernel
  house.svm <- svm(MEDV ~ ., data = trainset, kernel = 'radial', cost = 1000, gamma = 0.001)
  svm.pred <- predict(house.svm, testset[,-14])
  MSE.svm[i] <- crossprod(svm.pred - testset[,14]) / length(testindex)
  R2.svm[i] <- cor(svm.pred, testset[,14])^2
  
  # SVR using a linear kernel
  house.svml <- svm(MEDV ~ ., data = trainset, kernel = 'linear', cost = 7.4, epsilon = 0.1)
  svml.pred <- predict(house.svml, testset[,-14])
  MSE.svml[i] <- crossprod(svml.pred - testset[,14]) / length(testindex)
  R2.svml[i] <- cor(svml.pred, testset[,14])^2
  
  # mlr
  house.mlr <- lm(MEDV ~ ., data = trainset)
  mlr.pred <- predict(house.mlr, testset[,-14], interval = "prediction")
  MSE.mlr[i] <- crossprod(mlr.pred[,1] - testset[,14]) / length(testindex)
  R2.mlr[i] <- cor(mlr.pred[,1], testset[,14])^2
  
  # ridge
  house.ridge <- lm.ridge(MEDV ~ ., data = trainset, lambda = seq(0, 10, 0.01))
  best.newlambda = as.numeric(names(which.min(house.ridge$GCV)))
  house.ridge <- lm.ridge(MEDV ~ ., data = trainset, lambda = best.newlambda)
  ridge.pred <- coef(house.ridge)[1]
  for (j in 2:length(coef(house.ridge))){
    ridge.pred <- ridge.pred + coef(house.ridge)[j] * testset[,j-1]
  }
  MSE.ridge[i] <- crossprod(ridge.pred - testset[,14]) / length(testindex)
  R2.ridge[i] <- cor(ridge.pred, testset[,14])^2
  
  # lasso
  house.las <- cv.glmnet(x.train, trainset[, 14], alpha = 1)
  best.lam <- house.las$lambda.min
  alter.lam <- house.las$lambda.1se
  x.test <- data.matrix(testset[,-14])
  lasso.predmin <- predict(house.las, s = best.lam, newx = x.test)
  MSE.lassomin[i] <- crossprod(lasso.predmin - testset[,14]) / length(testindex)
  R2.lassomin[i] <- cor(lasso.predmin, testset[,14])^2
  lasso.pred <- predict(house.las, s = alter.lam, newx = x.test)
  MSE.lasso[i] <- crossprod(lasso.pred - testset[,14]) / length(testindex)
  R2.lasso[i] <- cor(lasso.pred, testset[,14])^2
}

boxplot(MSE.svm, MSE.svml, MSE.mlr, MSE.ridge, MSE.lassomin, MSE.lasso)
boxplot(R2.svm, R2.svml, R2.mlr, R2.ridge, R2.lassomin, R2.lasso)

jpeg('mse.jpg')
boxplot(MSE.mlr, MSE.ridge, MSE.lassomin, MSE.svml, MSE.svm, 
        names = c('MLR', 'ridge', 'lasso', 'SVR (lin)', 'SVR (RBF)'), 
        xlab = 'model', ylab = 'MSE')
dev.off()

jpeg('R2.jpg')
boxplot(R2.mlr, R2.ridge, R2.lassomin, R2.svml, R2.svm, 
        names = c('MLR', 'ridge', 'lasso', 'SVR (lin)', 'SVR (RBF)'),
        xlab = 'model', ylab = 'R^2')
dev.off()

mean(MSE.svm)
mean(MSE.svml)
mean(MSE.mlr)
mean(MSE.ridge)
mean(MSE.lassomin)
mean(R2.svm)
mean(R2.svml)
mean(R2.mlr)
mean(R2.ridge)
mean(R2.lassomin)

house.las <- cv.glmnet(x.train, trainset[, 14])
house.las$lambda.min

house.las <- cv.glmnet(x.train, trainset[, 14], alpha = 1)
best.lam <- house.las$lambda.min


x.train <- data.matrix(trainset[,-14])
house.lasso <- lars(x.train, trainset[,14], type = 'lasso')
house.newCp <- summary(house.lasso)$Cp
best.newCp <- (1:length(house.newCp))[house.newCp == min(house.newCp)]
coef.lasso = coef(house.lasso)[best.newCp,]
predict.lars(house.lasso, testset[, -14], s = best.newCp, type = 'fit')

house.ridge <- lm.ridge(MEDV ~ ., data = trainset, lambda = seq(0, 10, 0.01))
best.newlambda = as.numeric(names(which.min(house.ridge$GCV)))
best.newlambda
house.ridge <- lm.ridge(MEDV ~ ., data = trainset, lambda = best.newlambda)


tune.svm(testset[,-14], testset[,14], kernel = 'radial', cost = 10^(0:3), gamma = c(0.0001, 0.001, 0.01))

tune.svm(testset[,-14], testset[,14], kernel = 'linear', cost = 10^(0:3), epsilon = c(0.1, 1))
data(Ozone, package="mlbench")

## best parameters fixed to be c = 7.4 and epsilon = 0.1 (default)

## best parameters fixed to be c = 1000, gamma = 0.001 and epsilon = (default)
svm.house <- svm(MEDV ~ ., data = house.data, kernel = 'radial', cost = 1000, gamma = 0.001)
svm.house$coefs
length(house.data[,1])
w <- t(svm.house$coefs) %*% svm.house$SV 
b <- svm.house$rho 
house.train <- read.table('housing_train.txt', header = F)
colnames(house.train) <- c('CRIM', 'ZN', 'INDUS', 'CHAS', 'NOX', 'RM', 'AGE', 'DIS', 'RAD', 'TAX', 'PTRATIO', 'B', 'LSTAT', 'MEDV')
house.test <- read.table('housing_test.txt', header = F)
colnames(house.test) <- c('CRIM', 'ZN', 'INDUS', 'CHAS', 'NOX', 'RM', 'AGE', 'DIS', 'RAD', 'TAX', 'PTRATIO', 'B', 'LSTAT', 'MEDV')

# mlr using the training set
house.lm2 <- lm(MEDV ~ ., data = house.train)
layout(matrix(c(1,2,3,4),2,2))
plot(house.lm2)
summary(house.lm2)
# house.lm2 is equivalent to house.lm3
house.lm3 <- lm(MEDV ~ CRIM + ZN + INDUS + factor(CHAS) + NOX + RM + AGE + DIS 
                + RAD + TAX + PTRATIO + B + LSTAT, data = house.train)
summary(house.lm3)
house.predict1 <- predict(house.lm2, x.test, interval = "prediction")

MSE.predict1 <- sum((house.test[,14] - house.predict1[,1])^2) / length(house.test[,14])

# ?
house.SVR1 <- read.table('housing_predict1', header = F)
MSE.SVR1 <- sum((house.test[,14] - house.SVR1)^2) / length(house.test[,14])

str(house.data)
summary(house.data)

install.packages("e1071")
library('e1071')
# linear kernel
house.svm1 <- svm(MEDV ~ ., kernel = 'linear', data = house.train)
summary(house.svm1)
# RBF kernel
house.svm2 <- svm(MEDV ~ ., data = house.train)
summary(house.svm2)

# test data without response
x.test <- subset(house.test, select = -MEDV)

house.pred1 <- predict(house.svm1, x.test)
house.pred2 <- predict(house.svm2, x.test)
MSE.pred1 <- sum((house.pred1 - house.test$MEDV)^2) / length(house.pred1)
MSE.pred2 <- sum((house.pred2 - house.test$MEDV)^2) / length(house.pred2)

# 5-fold cross validation
house.svm3 <- svm(MEDV ~ ., kernel = 'linear', cross = 5, data = house.train)
summary(house.svm3)

# scale? MSE.pred4 = MSE.pred1
house.svm4 <- svm(MEDV ~ ., kernel = 'linear', scale = c(rep(TRUE, 13), FALSE), data = house.train)
summary(house.svm4)
house.pred4 <- predict(house.svm4, x.test)
MSE.pred4 <- sum((house.pred4 - house.test$MEDV)^2) / length(house.pred4)

# MSE.pred5 slightly greater than MSE.pred1
house.svm5 <- svm(MEDV ~ ., kernel = 'linear', scale = c(FALSE, rep(TRUE, 13)), data = house.train)
summary(house.svm5)
house.pred5 <- predict(house.svm5, x.test)
MSE.pred5 <- sum((house.pred5 - house.test$MEDV)^2) / length(house.pred5)

##use house.svm1
house.tune1 <- tune.svm(MEDV~., data = house.train, kernel = 'linear', cost = 2^(1:3), epsilon = 0.1^(-1:1))
house.tune2 <- tune.svm(MEDV~., data = house.train, kernel = 'linear', cost = 2^(1:3), epsilon = 0.1^(-2:0))
house.svm3 <- svm(MEDV ~ ., kernel = 'linear', cross = 10, cost = 4, data = house.train)
summary(house.svm3)

#best parameters: cost = 7.2 and epsilon = 0.275
tune.svm(MEDV~., data = house.train, kernel = 'linear', cost = seq(1, 8, by = 0.2), 
         epsilon = seq(0.1, 0.4, by = 0.025))
house.svm6 <- svm(MEDV ~ ., kernel = 'linear', cost = 7.2, epsilon = 0.275, data = house.train)
summary(house.svm6)
house.pred6 <- predict(house.svm6, x.test)
MSE.pred6 <- sum((house.pred6 - house.test$MEDV)^2) / length(house.pred6)

r2.pred6 <- sum((house.pred6 - mean(house.test$MEDV))^2) / sum((house.test$MEDV - mean(house.test$MEDV))^2)
r2.pred1 <- sum((house.pred1 - mean(house.test$MEDV))^2) / sum((house.test$MEDV - mean(house.test$MEDV))^2)
r2.predict1 <- sum((house.predict1[,1] - mean(house.test$MEDV))^2) / sum((house.test$MEDV - mean(house.test$MEDV))^2)
summary(house.lm2)$r.squared
cor(fitted(house.lm2), house.train$MEDV)^2
cor(house.pred6, house.test$MEDV)^2
cor(house.pred1, house.test$MEDV)^2
cor(house.predict1[,1], house.test$MEDV)^2
cor(house.pred2, house.test$MEDV)^2
c <- max(house.train$MEDV) - min(house.train$MEDV)
house.tune2 <- tune.svm(MEDV~., data = house.train, kernel = 'radial', gamma = 2.^(1:7), epsilon = 1:5, cross = 5)
house.tune3 <- tune(svm, MEDV~., data = house.train, type = 'eps-regression'
            ranges = list(gamma = 2^(-7:7), cost = c, epsilon = 0:5),
            tunecontrol = tune.control(sampling = "fix")
)
x.train = subset(house.train, select = -MEDV)
y.train = house.train$MEDV
house.tune2 <- tune.svm(x.train, y.train, kernel = 'radial', gamma = 2.^(1:7), epsilon = 1:5, cross = 5)

r2.lm2 <- sum((fitted(house.lm2) - mean(house.train$MEDV))^2) / sum((house.train$MEDV - mean(house.train$MEDV))^2)

## Bayesian linear model
library(BMS)
house.zlm1 <- zlm(MEDV ~ ., data = house.train, g = "UIP")
house.zpred1 <- predict(house.zlm1, x.test, interval = "prediction")
MSE.zpred1 <- sum((house.test[,14] - house.zpred1)^2) / length(house.test[,14])
cor(house.zpred1, house.test$MEDV)^2
cor(house.zpred1, house.test[,14])^2

jpeg('zlm_pred')
plot(house.test$MEDV, house.zpred1, ann = F, pch = 19, col = 'blue')
title(main = 'Test set', xlab = 'observed price', ylab = 'predicted price')
lines(x,y=x, lwd = 2)
dev.off()

jpeg('zlm_train')
plot(house.train$MEDV, fitted(house.zlm1), ann = F, pch = 19, col = 'blue')
title(main = 'Training set', xlab = 'observed price', ylab = 'predicted price')
lines(x,y=x, lwd = 2)
dev.off()

jpeg('mlr_pred')
plot(house.test$MEDV, house.predict1[,1], ann = F, pch = 19, col = 'blue')
title(main = 'Test set', xlab = 'observed price', ylab = 'predicted price')
lines(x,y=x, lwd = 2)
dev.off()

jpeg('mlr_train')
plot(house.train$MEDV, fitted(house.lm2), ann = F, pch = 19, col = 'blue')
title(main = 'Training set', xlab = 'observed price', ylab = 'predicted price')
lines(x,y=x, lwd = 2)
dev.off()

jpeg('svr_pred')
plot(house.test$MEDV, house.pred1, ann = F, pch = 19, col = 'blue')
title(main = 'Test set', xlab = 'observed price', ylab = 'predicted price')
lines(x,y=x, lwd = 2)
dev.off()

jpeg('svr_train')
plot(house.train$MEDV, fitted(house.svm1), ann = F, pch = 19, col = 'blue')
title(main = 'Training set', xlab = 'observed price', ylab = 'predicted price')
lines(x,y=x, lwd = 2)
dev.off()

# correlation in the data
library(lattice)
levelplot(cor(house.data[,-14]))
cor(house.data[,-14])
