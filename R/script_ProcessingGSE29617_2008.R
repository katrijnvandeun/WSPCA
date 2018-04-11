
setwd('C:/Users/Katrijn/OneDrive/Documents/Manuscripts/InProgress/WSPCA/R')


############################################
#             NEEDED PACKAGES              #
############################################

packages <- c("limma", "statmod", "Biobase", "affy", "affyPLM","GEOquery","LMGene")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(setdiff(packages, rownames(installed.packages())))  
}
 
library(affy)
library(statmod)
library(limma)
library(Biobase)
library(GEOquery)
library(affyPLM)
install.packages("GESTr") 
library(GESTr)
library(LMGene)

accnr<-c("GSE29617")#enter here the accession number

# load series and platform data from GEO
getGEOSuppFiles(accnr)
setwd(paste("./",accnr,sep=""))
untar(paste(accnr,"_RAW.tar",sep=""), exdir="rawdata")
setwd("./rawdata")
cels <- list.celfiles()
data<-ReadAffy(filenames=cels)
#show(data): 80 samples; annotation=hthgu133pluspm
eset <- rma(data)
eset<-(exprs(eset))

#estimate additive and multiplicative error variance

#fitRockeDurbin(x,theta)#x:expression set; theta: does not influence pars
RLpars<-fitRockeDurbin(eset,0)
RLparsD0<-fitRockeDurbin(eset[,samplecov$characteristics_ch1.2=="time point: D0"],0)
RLparsD3<-fitRockeDurbin(eset[,samplecov$characteristics_ch1.2=="time point: D3"],0)

#muest<-exprs(eset)#-RLpars$alpha
weights<-(eset^2)/(eset^2*RLpars$sd_epsilon^2+RLpars$sd_epsilon^2)
plot(eset,weights)#takes quite a lot of time

write.table(eset,sep = "\t", file=(paste(accnr,"_rma.txt",sep="")),quote=FALSE,row.names = FALSE,col.names = FALSE)
write.table(RLparsD0,sep = "\t", file=(paste(accnr,"_RLparsD0.txt",sep="")),quote=FALSE,row.names = FALSE,col.names = FALSE)
write.table(RLparsD3,sep = "\t", file=(paste(accnr,"_RLparsD3.txt",sep="")),quote=FALSE,row.names = FALSE,col.names = FALSE)

#sample information
gset <- getGEO(paste(accnr), GSEMatrix = TRUE)
if (length(gset) > 1) idx <- grep(paste(accnr), attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
colnames(pData(gset))
head(pData(gset))
samplecov<-pData(gset)
colnames(samplecov)
forexport<-samplecov[,c("characteristics_ch1.1","characteristics_ch1.2")]
forexport_titers<-samplecov[,c("characteristics_ch1.4","characteristics_ch1.5","characteristics_ch1.6","characteristics_ch1.7",
"characteristics_ch1.8","characteristics_ch1.9")]
 

write.table(forexport,sep = "\t", file=(paste(accnr,"_pdata.txt",sep="")),quote=FALSE,row.names = FALSE,col.names = FALSE)
write.table(forexport_titers,sep = "\t", file=(paste(accnr,"_titers.txt",sep="")),quote=FALSE,row.names = FALSE,col.names = FALSE)