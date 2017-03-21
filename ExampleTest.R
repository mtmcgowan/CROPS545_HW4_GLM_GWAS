# 

mdpN_raw=read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
mdpN <- mdp_numeric[,2:3094]
row.names(mdpN) <- mdpN_raw[,1]

mdpPheno_raw=read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)
mdpCov_raw=read.table(file='http://zzlab.net/GAPIT/data/CROP545_Covariates.txt', head=T)