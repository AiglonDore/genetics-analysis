library(gglasso)
load("project_park.RData")

df <- as.data.frame(Xmat)
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
sortedBetas <- sort(x = absBetas, index.return = TRUE, decreasing = TRUE)