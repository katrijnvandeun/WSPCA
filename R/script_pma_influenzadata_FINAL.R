rm(list = ls())

####Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#in R: setwd(getSrcDirectory()[1])
#setwd('C:/Users/kvandeun/surfdrive/Manuscripts/InProgress/WSPCA/R')
#####

############################################
#             NEEDED PACKAGES              #
############################################

packages <- c("Biobase", "annotate", "hgu133plus2.db", "impute")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(setdiff(packages, rownames(installed.packages())))  
}
packages <- c("PMA",'MASS')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library(hgu133plus2.db)
library(impute)
library(PMA)
library(MASS)

#Get official gene symbols: TAKES SOME TIME!
x <- hgu133plus2SYMBOL 
symbolall<-Lkeys(x)
for (i in 1:length(symbolall))
{
  symbolall[i] <- get(symbolall[i],env=hgu133plus2SYMBOL)
}


#####################################
#
#      Extract data
#
#####################################

#load 'training' sample
DATA2008<-read.table(paste('../DATA/DATA2008.txt',sep=""))
DATA2008<-t(DATA2008)
DATA2008 <- DATA2008[,1:54675]
TITER2008<-read.table(paste('../DATA/TIVtiter2008.txt',sep=""))
#load independent test sample, 2007 season
# Another micro-array platform was used in that season, no control probes
DATA2007<-read.table(paste('../DATA/DATA2007.txt',sep=""))
DATA2007<-t(DATA2007)
TITER2007<-read.table(paste('../DATA/TIVtiter2007.txt',sep=""))
#load 'dependent' variable for both seasons and transform (see Nakaya et al., 2011)
y<-log2(TITER2008[,4])
y2007<-log2(TITER2007[,4])

#####################################
#
#      PMA analysis
#
#####################################

#set parameters for the PMA analysis
R<-3
L_627 = (sqrt(101)) #upperbound on number of non-zero coefficients is 627 over three components
result<-SPC(DATA2008, sumabsv=L_627, niter=50, K=R, orth=TRUE, trace=TRUE, v=NULL,
            center=FALSE, cnames=NULL, vpos=FALSE, vneg=FALSE, compute.pve=TRUE)
sum(result$v!=0)
colSums(result$v!=0)

#calculation of fit
T<-result$u#orthonormal component scores as in WSPCA
df<-data.frame(cbind(y,T))
colnames(df)<-c("y","t1","t2",'t3')
regr<-lm(y~t1+t2+t3,data=df)
rsq2008<-summary(regr)$r.squared
rsq2008
MSEtrain <- sum((y-regr$fitted.values)^2)/length(y)
MSEtrain

#calculation of prediction error
Pt <- diag(result$d, nrow=R)%*%t(result$v)
T2007<-DATA2007%*%ginv(Pt)#From X=USV'+E in PMA it follows that U=X(SV')^+
newdf<-as.data.frame(cbind(y2007,T2007))
colnames(newdf)<-c("y","t1","t2","t3")
regr2007<-predict(regr,newdf,type="response")
cor(y2007,regr2007)^2
MSEtest <- sum((y2007-regr2007)^2)/length(y2007)
MSEtest

#####################################
#
#      WRITE output
#
#####################################


#store official gene symbols and associated loading for selected probe sets
for (r in 1:R) {
  df = cbind(symbolall[result$v[,r]!=0],result$v[result$v[,r]!=0,r])
  write.table(df, paste("../RESULTS/GENEIDS_pma_R",r,"_.txt",sep = ""), row.names=F, col.names=F,sep=" ")
}

#for the WSPCA results obtained with our own procedure in MATLAB
WSPCAloadings<-read.table(paste('../DATA/wspcaloadings.txt',sep=""))
sum(WSPCAloadings!=0)
DF <- c()
for (r in 1:R) {
  df = cbind(symbolall[WSPCAloadings[,r]!=0],WSPCAloadings[WSPCAloadings[,r]!=0,r])
  DF = rbind(DF,df)
}
write.table(DF, paste("../RESULTS/SelectedProbesets/GENEIDS_wspca_R",'ALL',"_.txt",sep = ""), row.names=F, col.names=F,sep=" ")

