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

variables <- df[sortedIndex[1:10000]]

variables$Y <- df$Y
variables <- na.omit(variables)

tvariables <- t(variables)
kmeans.res <- kmeans(x = tvariables[-c(10001), ], centers = 50, iter.max = 100)
gglasso.res <- gglasso(as.matrix(variables[-c(10001)]), variables[, 10001], group = kmeans.res$cluster)
plot(gglasso.res)

library(glmnet)
lambda <- cv.gglasso(as.matrix(variables[-c(10001)]), variables[, 10001], group = kmeans.res$cluster, nfolds = 10)
lambda.min <- lambda$lambda.min

gglasso.opt <- gglasso(as.matrix(variables[-c(10001)]), variables[, 10001], group = kmeans.res$cluster, lambda = lambda.min)