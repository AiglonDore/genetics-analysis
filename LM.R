# Loading data
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

x <- matrix(data = kmeans(Xmat, 50), nrow = 413)
df <- as.data.frame(x)
rm(x)
for (i in 1:ncol(df)){
    df[, i] <- as.numeric(df[, i])
}
rm(i)
df$Y <- pheno.df$Protein.content
library(flexmix)
results <- list(c(),c())
names(results) <- c("r2","bic")
for (i in 1:20){
    reg <- lm(Y~.,data = df)
    #Variable selection
    df <- na.omit(df)
    vselect <- step(reg, direction = "both")
    results$r2 <- append(results$r2, summary(vselect)$r.squared)
    results$bic <- append(results$bic, BIC(vselect))
}

print(paste("Moyenne des R2: ",mean(results$r2)))
print(paste("Maximum des R2: ",max(results$r2)))
print(paste("Moyenne des BIC: ",mean(results$bic)))
print(paste("Minimum des BIC: ",max(results$bic)))

#K-fold
library(caret)
kfold <- trainControl(method = "cv", number = 25)
model <- train(Y~., data = df, method = "lm", trControl = kfold)
summary(model)
