#sparse PCA with missing data using PMA package of Witten et al.

setwd('C:/Users/kvandeun/surfdrive/Manuscripts/InProgress/WSPCA/R')
source("https://bioconductor.org/biocLite.R")
biocLite("impute")
library(impute)
library(PMA)
install.packages("psych", type="binary")
library(psych)

#######ONLY RUN IF REANALYSIS NEEDED
RESULTPrandom <- c()
for (teller in (1:720))
{
	DATA<-read.table(paste('../DATA/DATAnomissb',teller,'.dat',sep=""))
	TRUEDATA<-read.table(paste('../DATA/TRUEnomissb',teller,'.dat',sep=""))
	TRUEd<-as.matrix(TRUEDATA)
	W<-read.table(paste('../DATA/WEIGHTSnomissb',teller,'.dat',sep=""))
	DATAmiss<-as.matrix(DATA)
	DATAmiss[W==0]<-NaN
	DATAmiss[is.infinite(DATAmiss)]<-NaN
	PTRUE<-read.table(paste('../DATA/TRUEPnomissb',teller,'.dat',sep=""))
	PTRUE<-as.matrix(PTRUE)
	TTRUE<-read.table(paste('../DATA/TRUETnomissb',teller,'.dat',sep=""))
	TTRUE<-as.matrix(TTRUE)
	R<-dim(PTRUE)[2]
	J<-dim(PTRUE)[1]
	d<-dim(W)
	nrnonzeros_true <- J #both sparse and non-sparse generated data analyzed assuming half of the R*J loadings are zero
	result<-SPC(DATAmiss, sumabsv=sqrt(nrnonzeros_true/R), niter=50, K=R, orth=TRUE, trace=TRUE, v=NULL,
	center=FALSE, cnames=NULL, vpos=FALSE, vneg=FALSE, compute.pve=TRUE)
	nrzeros_spc<-sum(rowSums(result$v==0))
	Wmiss<-W
	Wmiss[W!=0]<-1
	Wmiss[is.infinite(DATAmiss)]<-0
	DATAhat<-result$u%*%diag(result$d,nrow=R)%*%t(result$v)
	devtrue<-(Wmiss)*(TRUEd-DATAhat)
	WTRUE<-(Wmiss)*TRUEd
	presstrueiw<-sum(colSums(devtrue^2))/sum(colSums(WTRUE^2))
	congrmat=psych::factor.congruence(PTRUE,result$v%*%diag(result$d,nrow=R))
	spcr_tucongrP=max(abs(congrmat[1,1])+abs(congrmat[2,2]),abs(congrmat[2,1])+abs(congrmat[1,2]))/2
	congrmat=psych::factor.congruence(TTRUE,result$u)
	spcr_tucongrT=max(abs(congrmat[1,1])+abs(congrmat[2,2]),abs(congrmat[2,1])+abs(congrmat[1,2]))/2
	RESULTPrandom<-rbind(RESULTPrandom,c(presstrueiw,spcr_tucongrT,spcr_tucongrP))
}
write.table(RESULTPrandom,file="../RESULTS/PMDresultPrandom.txt",sep="\t",eol="\n")
#######


#RESULTPrandom<-read.table("../RESULTS/PMDresultPrandom.txt",sep="\t")
RESULT <- RESULTPrandom
######BOX PLOTS ##############


#load PCA,wPCA,wSPCA results
matlabresults=read.table("RESULTnomiss_sqrtW_Prandom.txt",sep="\t")
colnames(matlabresults) <- c('ADD','MULT','PrMiss','PrSp','I','J','R',rep('Press',3),rep('Tu_T',3),'Tu_Pca','Tu_WPCA','Tu_WSPCA',
                             'W0','W20','W50','W90')

#define method factor
method=rep(c(" pca"," wpca","spca","wspca"),c(360,360,360,360))

#SPARSE
nonsparse <- which(matlabresults[,4]==0.5)
# define factors of the simulation experiment
#noiseadd noisemult propmiss propsparse size(DATA,1) size(DATA,2) 
adderr<-as.factor(rep(matlabresults[nonsparse,1],4))
multerr<-as.factor(rep(matlabresults[nonsparse,2],4))
I<-as.factor(rep(matlabresults[nonsparse,5],4))
J<-as.factor(rep(matlabresults[nonsparse,6],4))

#define performance measures
tuckerP=c(matlabresults[nonsparse,14],matlabresults[nonsparse,15],RESULT[nonsparse,3],matlabresults[nonsparse,16])
tuckerT=c(matlabresults[nonsparse,11],matlabresults[nonsparse,12],RESULT[nonsparse,2],matlabresults[nonsparse,13])


#TUCKER CONGRUENCE LOADINGS
Level_I <-'100'
Level_J <-'1000'
  
pdf(paste("boxplot_Prandom_NONsparse_congruenceP_I",Level_I,"_J",Level_J,".pdf",sep = ""),width=12,height = 8,paper="special")
par(mfrow=c(1,3))

a=boxplot(tuckerP[multerr=="0" & J==Level_J & I==Level_I]~method[multerr=="0"& J==Level_J & I==Level_I]+adderr[multerr=="0"& J==Level_J & I==Level_I],
          at=c(0.1,1.3,2.5,3.7,    
               7.5,8.7,9.9,11.1,    
               14.9,16.1,17.3,18.5),
          xaxt="n",ylab="Tucker Congruence: loadings",main="MULT=0",ylim=c(0,1),col=c(1,2,3,4))
axis(1, at = c(1.8 , 9.3 , 16.7), labels = c("ADD=0","ADD=0.01","ADD=0.05") , tick=FALSE , cex=0.3)
abline(v=c(5.6,13),lty=1, col="grey")
# Add a legend
legend("bottomright", legend = c("PCA", "WPCA", "SPCA","WSPCA"), col=c(1 ,2,3,4),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = 0.01)

a=boxplot(tuckerP[multerr=="0.01" & J==Level_J & I==Level_I]~method[multerr=="0.01"& J==Level_J & I==Level_I]+adderr[multerr=="0.01"& J==Level_J & I==Level_I],
          at=c(0.1,1.3,2.5,3.7,    
               7.5,8.7,9.9,11.1,    
               14.9,16.1,17.3,18.5),
          xaxt="n",ylab="Tucker Congruence: loadings",main="MULT=0.01",ylim=c(0,1),col=c(1,2,3,4))
axis(1, at = c(1.8 , 9.3 , 16.7), labels = c("ADD=0","ADD=0.01","ADD=0.05") , tick=FALSE , cex=0.3)
abline(v=c(5.6,13),lty=1, col="grey")
# Add a legend
legend("bottomright", legend = c("PCA", "WPCA", "SPCA","WSPCA"), col=c(1 ,2,3,4),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = 0.01)

a=boxplot(tuckerP[multerr=="0.05" & J==Level_J & I==Level_I]~method[multerr=="0.05"& J==Level_J & I==Level_I]+adderr[multerr=="0.05"& J==Level_J & I==Level_I],
          at=c(0.1,1.3,2.5,3.7,    
               7.5,8.7,9.9,11.1,    
               14.9,16.1,17.3,18.5),
          xaxt="n",ylab="Tucker Congruence: loadings",main="MULT=0.05",ylim=c(0,1),col=c(1,2,3,4))
axis(1, at = c(1.8 , 9.3 , 16.7), labels = c("ADD=0","ADD=0.01","ADD=0.05") , tick=FALSE , cex=0.3)
abline(v=c(5.6,13),lty=1, col="grey")
# Add a legend
legend("bottomright", legend = c("PCA", "WPCA", "SPCA","WSPCA"), col=c(1 ,2,3,4),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = 0.01)

dev.off()



#PRESS
Level_I <-'100'
Level_J <-'1000'
press=c(matlabresults[nonsparse,8],matlabresults[nonsparse,9],RESULT[nonsparse,1],matlabresults[nonsparse,10])
press <- log(press+1)
m <- min(press)
M <- 0.00001#2# max(press)# 2 #outliers present 

pdf(paste("boxplot_Prandom_NONsparse_press_I",Level_I,"_J",Level_J,".pdf",sep = ""),width=12,height = 8,paper="special")

par(mfrow=c(1,3))

a=boxplot(press[multerr=="0" & J==Level_J & I==Level_I]~method[multerr=="0"& J==Level_J & I==Level_I]+adderr[multerr=="0"& J==Level_J & I==Level_I],
          at=c(0.1,1.3,2.5,3.7,    
               7.5,8.7,9.9,11.1,    
               14.9,16.1,17.3,18.5),
          xaxt="n",ylab="Badness-of-fit",main="MULT=0",ylim=c(m,M),col=c(1,2,3,4))
axis(1, at = c(1.8 , 9.3 , 16.7), labels = c("ADD=0","ADD=0.01","ADD=0.05") , tick=FALSE , cex=0.3)
abline(v=c(5.6,13),lty=1, col="grey")
# Add a legend
legend("topleft", legend = c("PCA", "WPCA", "SPCA","WSPCA"), col=c(1 ,2,3,4),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = 0.01)

a=boxplot(press[multerr=="0.01" & J==Level_J & I==Level_I]~method[multerr=="0.01"& J==Level_J & I==Level_I]+adderr[multerr=="0.01"& J==Level_J & I==Level_I],
          at=c(0.1,1.3,2.5,3.7,    
               7.5,8.7,9.9,11.1,    
               14.9,16.1,17.3,18.5),
          xaxt="n",ylab="Badness-of-fit",main="MULT=0.01",ylim=c(m,M),col=c(1,2,3,4))
axis(1, at = c(1.8 , 9.3 , 16.7), labels = c("ADD=0","ADD=0.01","ADD=0.05") , tick=FALSE , cex=0.3)
abline(v=c(5.6,13),lty=1, col="grey")
# Add a legend
legend("topleft", legend = c("PCA", "WPCA", "SPCA","WSPCA"), col=c(1 ,2,3,4),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = 0.01)

a=boxplot(press[multerr=="0.05" & J==Level_J & I==Level_I]~method[multerr=="0.05"& J==Level_J & I==Level_I]+adderr[multerr=="0.05"& J==Level_J & I==Level_I],
          at=c(0.1,1.3,2.5,3.7,    
               7.5,8.7,9.9,11.1,    
               14.9,16.1,17.3,18.5),
          xaxt="n",ylab="Badness-of-fit",main="MULT=0.05",ylim=c(m,M),col=c(1,2,3,4))
axis(1, at = c(1.8 , 9.3 , 16.7), labels = c("ADD=0","ADD=0.01","ADD=0.05") , tick=FALSE , cex=0.3)
abline(v=c(5.6,13),lty=1, col="grey")
# Add a legend
legend("topleft", legend = c("PCA", "WPCA", "SPCA","WSPCA"), col=c(1 ,2,3,4),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = 0.01)

dev.off()

