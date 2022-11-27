#Loading data
load("project_park.RData")

#Building dataframe
df <- as.data.frame(Xmat)
df$Y <- pheno.df$Protein.content

#Removing useless data
rm(geno.df)
rm(pheno.df)
rm(Xmat)
gc()