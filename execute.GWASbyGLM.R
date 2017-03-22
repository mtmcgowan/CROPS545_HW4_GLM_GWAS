mdpN_raw=read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
mdpN <- mdpN_raw[,2:3094]
row.names(mdpN) <- mdpN_raw[,1]

mdpPheno_raw<-read.table(file="http://zzlab.net/GAPIT/data/CROP545_Phenotype.txt",head=T)
mdpCov_raw<-read.table(file='http://zzlab.net/GAPIT/data/CROP545_Covariates.txt', head=T)

myGM<-read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)

source("GWASbyGLM2.R")
source("GLM.R")
source("G2P.R")
source("GWASbyCor.R")

geno_file<-mdpN[,-1]
index1to5<-myGM[,2]<6
X1to5<-mdpN[,index1to5]

set.seed(99164)
mySim<-G2P(X= X1to5,h2=.75,alpha=1,NQTN=10,distribution="norm")#phenotype simulation

###Question#4
myGLM<-GWASbyGLM2(mdpN=mdpN,mdpPheno=mdpPheno_raw[,-1],mdpCov=mdpCov_raw[,-1],option="cov")


###QQ plot
p.obs=myGLM
m2=length(p.obs)
p.uni=runif(m2,0,1)
order.obs=order(p.obs)
order.uni=order(p.uni)
plot(-log10(p.uni[order.uni]), -log10(p.obs[order.obs]), 
     xlab="Expected -log10(p)", ylab="Observed -log10(p)")
abline(a=0, b=1, col="red")

##Manhattan plot##

#bonferroni cutoff
cutoff=0.05
P.value=myGM
order.SNP=order(P.value)
bonfcutoff <- cutoff/m

color.vector <- rep(c("deepskyblue","orange","forestgreen","indianred3"),10)
m=nrow(myGM)
plot((-log10(myGLM))~seq(1:m),col=color.vector[myGM[,2]])
abline(h=-log10(bonfcutoff), col="red")


##Question#2
myGLM_PCA<-GWASbyGLM2(mdpN=mdpN,mdpPheno=mdpPheno_raw[,-1],mdpCov=mdpCov_raw[,-1],option="pca")

cutoff=0.05
P.value=myGM_PCA
order.SNP=order(P.value)
bonfcutoff <- cutoff/m

color.vector <- rep(c("deepskyblue","orange","forestgreen","indianred3"),10)
m=nrow(myGM)
plot((-log10(myGLM_PCA))~seq(1:m),col=color.vector[myGM[,2]])
abline(h=-log10(bonfcutoff), col="red")

##Question#5
###comparison between GWASbyCor and GWASbyGLM

source("GWASbyCor.R")
geno_file<-mdpN[,-1]
index1to5<-myGM[,2]<6
X1to5<-geno_file[,index1to5]

set.seed(99164)

mySim<-G2P(X= X1to5,h2=0.75,alpha=1,NQTN=10,distribution="norm")#phenotype simulation

myCor<-GWASbyCor(X=mdpN,y=mySim$y)
myGLM_PCA<-GWASbyGLM2(mdpN=mdpN,mdpPheno=mySim$y,mdpCov=mdpCov_raw[,-1],option="pca")

#QQ plot

QQ_plot_cor=function(p.obs){
  m2=length(p.obs)
  p.uni=runif(m2,0,1)
  order.obs=order(p.obs)
  order.uni=order(p.uni)
  plot(-log10(p.uni[order.uni]), -log10(p.obs[order.obs]), 
       xlab="Expected", ylab="Observed",main="QQ plot for GWASbyCor")
  abline(a=0, b=1, col="red")
}
QQ_plot_glm=function(p.obs){
  m2=length(p.obs)
  p.uni=runif(m2,0,1)
  order.obs=order(p.obs)
  order.uni=order(p.uni)
  plot(-log10(p.uni[order.uni]), -log10(p.obs[order.obs]), 
       xlab="Expected", ylab="Observed",main="QQ plot for GWASbyGLM")
  abline(a=0, b=1, col="red")
}
QQplot_GWASbyCor=QQ_plot_cor(p.obs=myCor)
QQplot_GWASbyGLM=QQ_plot_glm(p.obs=myGLM_PCA)



##Manhattan plot##

Manhattan_plot_cor=function(p){
color.vector <- rep(c("deepskyblue","orange","forestgreen","indianred3"),10)
m=nrow(myGM)
plot(t(-log10(p))~seq(1:m),col=color.vector[myGM[,2]],main="Plot for GWASbyCor")
abline(v=mySim$QTN.position, lty = 2, lwd=2, col = "black")
}


Manhattan_plot_glm=function(p){
color.vector <- rep(c("deepskyblue","orange","forestgreen","indianred3"),10)
m=nrow(myGM)
plot((-log10(p))~seq(1:m),col=color.vector[myGM[,2]],main="Plot for GWASbyGLM")
abline(v=mySim$QTN.position, lty = 2, lwd=2, col = "black")
}
Manhattan_cor=Manhattan_plot_cor(p=myCor)##GWASbyCor
Manhattan_cor=Manhattan_plot_glm(p=myGLM_PCA)###GWASbyGLM


###calculate the FDR, power###
install.packages("compiler")
library(compiler) #required for cmpfun
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

##repeat the function, but with the same h2=0.75
statRep_cor=function(nrep,GM){
  
  stat_cor=replicate(nrep, {
    mySim=G2P(X=X1to5,h2=0.75,alpha=1,NQTN=10,distribution="norm")
    myCor= GWASbyCor(X=mdpN,y=mySim$y)
    seqQTN=mySim$QTN.position
    myGWAS=cbind(myGM,t(myCor),NA)
    myStat_cor=GAPIT.FDR.TypeI(WS=c(1e0,1e3,1e4,1e5), GM=myGM,seqQTN=mySim$QTN.position,GWAS=myGWAS,maxOut=100,MaxBP=1e10)
    
  })}

statRep_glm=function(nrep,GM){
  stat_glm=replicate(nrep, {
    mySim=G2P(X=X1to5,h2=0.75,alpha=1,NQTN=10,distribution="norm")
    myGLM_PCA<-GWASbyGLM2(mdpN=mdpN,mdpPheno=mySim$y,mdpCov=mdpCov_raw[,-1],option="pca")
    seqQTN=mySim$QTN.position
    myGWAS=cbind(myGM,myGLM_PCA,NA)
    myStat_PCA=GAPIT.FDR.TypeI(WS=c(1e0,1e3,1e4,1e5), GM=myGM,seqQTN=mySim$QTN.position,GWAS=myGWAS,maxOut=100,MaxBP=1e10)
  })
}

statRep_cor5=statRep_cor(nrep=5,GM=myGM) #repeat 5 times
statRep_glm5=statRep_glm(nrep=5,GM=myGM) #repeat 5 times,10 times take a long time

###GLM vs COR in terms of FDR and Power, Type1 error and power
par(mfrow=c(1,2),mar = c(5,2,5,2))
plot(statRep_glm5[[3]][,4],statRep_glm5[[2]],type="b",col="red",xlab="FDR",ylab="Power")
lines(statRep_cor5[[3]][,4],statRep_cor5[[2]],type="b",col="blue")
legend("topleft", c("GLM", "COR"), col=c("red", "blue"), lty=1, pch=c(1, 2), cex=0.8, bty="n")

plot(statRep_glm5[[4]][,4],statRep_glm5[[2]],type="b",col="red",xlab="Type1 error",ylab="Power")
lines(statRep_cor5[[4]][,4],statRep_cor5[[2]],type="b",col="blue")
legend("topleft", c("GLM", "COR"), col=c("red", "blue"), lty=1, pch=c(1, 2), cex=0.8, bty="n")


