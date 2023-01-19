load("project_park.RData")

kmeans_init <- function(d, k) {
    output <- list(floor(runif(n = k, min = 0, max = (d - 1))))
    return(output)
}

kmeans_assign <- function(data, centroids) {
    output <- list()
    for (i in 1:(length(centroids))) {
        x <- data[, i]
        dist <- list()
        for (j in 1:(length(centroids))) {
            y <- centroids[[j]]
            dist <- append(dist, sqrt(sum((x - y)^2)))
        }
        output <- append(output, which.min(dist))
    }
    return(output)
}

kmeans_update <- function(data, centroids, assignments) {
    output <- list()
    for (i in 1:(length(centroids))) {
        x <- data[, assignments == i]
        output <- append(output, colMeans(x))
    }
    return(output)
}

compare_lists <- function(a, b) {
    if (length(a) != length(b)) {
        return(FALSE)
    }
    for (i in 1:length(a)) {
        if (a[[i]] != b[[i]]) {
            return(FALSE)
        }
    }
    return(TRUE)
}

kmeans <- function(data, k, max_iter = NULL) {
    d <- ncol(data)
    index <- kmeans_init(d, k)
    centroids <- list()
    for (i in 1:(length(index))) {
        centroids <- append(centroids, data[, index[[i]]])
    }
    centroids2 <- list()
    i <- 0
    while (compare_lists(centroids, centroids2) || (!is.null(max_iter) && i < max_iter)) {
        i <- i + 1
        assignments <- kmeans_assign(data, centroids)
        centroids2 <- centroids
        centroids <- kmeans_update(data, centroids, assignments)
    }
    return(centroids)
}
r2 <- c()
bic <-c()
for (p in 1:20){
    x <- matrix(data = kmeans(Xmat, 50), nrow = 413)
    df <- as.data.frame(x)
    
    rm(x)
    for (i in 1:ncol(df)){
        df[, i] <- as.numeric(df[, i])
    }
    rm(i)
    df$Y <- pheno.df$Protein.content
    library(flexmix)
    names(results) <- c("r2","bic")
    reg <- lm(Y~.,data = df)
    #Variable selection
    df <- na.omit(df)
    vselect <- step(reg, direction = "both", trace = 0)
    r2 <- append(r2, summary(vselect)$r.squared)
    bicTest <- BIC(vselect)
    bic <- append(bic, bicTest)
    if (p == 1){
        bestBic <- bicTest
    }
    if (bicTest < bestBic){
        dfKfold <- df
        bestBic <- bicTest
    }
}

print(paste("Moyenne des R2: ",mean(r2)))
print(paste("Maximum des R2: ",max(r2)))
print(paste("Moyenne des BIC: ",mean(bic)))
print(paste("Minimum des BIC: ",min(bic)))

#K-fold
library(caret)
kfold <- trainControl(method = "cv", number = 25)
model <- train(Y~., data = dfKfold, method = "lm", trControl = kfold)
summary(model)

library(gglasso)
load("project_park.RData")

df <- as.data.frame(Xmat)
for (x in 1:nrow(df)) {
    for (y in 1:ncol(df)){
        if (is.na(df[x,y])) {
            df[x,y] = -1
        }
    }
}
df$Y <- pheno.df$Protein.content
rm(Xmat, geno.df, pheno.df)
gc()

betas = c()

for (i in 1:(ncol(df) - 1))
{
    reg <- lm(as.formula(paste("Y~",names(df)[i])), data = df)
    betas <- c(betas, reg$coefficients[2], recursive = TRUE)
}
rm(i, reg)
gc()

absBetas <- abs(betas)
sortedIndex <- order(x = absBetas, decreasing = TRUE)
sortedBetas <- c()
for (x in sortedIndex){
    sortedBetas <- append(sortedBetas, betas[x])
}
rm(x, absBetas)

library(dplyr)
library(glmnet)
library(flexmix)

variables50 <- df[sortedIndex[1:50]]
variables50$Y <- df$Y
variables50 <- na.omit(variables50)
lasso.cv <- cv.glmnet(as.matrix(variables50[-c(51)]), variables50[, 51])
plot(lasso.cv)
lasso.res <- glmnet(as.matrix(variables50[-c(51)]), variables50[, 51], lambda = lasso.cv$lambda)
plot(lasso.res)

lasso.res.coef <- predict.glmnet(lasso.res, type = "coefficients")
plot(lasso.res.coef)


tvariables50 <- t(variables50)
kmeans.res <- kmeans(x = tvariables50[-c(51), ], centers = 10, iter.max = 100)
gglasso.res <- gglasso(as.matrix(variables50[-c(51)]), variables50[, 51], group = kmeans.res$cluster)
plot(gglasso.res)
lambda <- cv.gglasso(as.matrix(variables50[-c(51)]), variables50[, 51], group = kmeans.res$cluster, nfolds = 5)
plot(lambda, main = "Choices of lambda")

gglasso.res.coef <- predict(gglasso.res, type = "class", newx = as.matrix(variables50[-c(51)]), mode = "lambda")
plot(gglasso.res.coef)