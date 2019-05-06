#########################################################################################
# A script to perform geneset enrichment analyis of genes identified by sPCA and wsPCA in 
# D3 data from Nakaya et al.
##########################################################################################

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../RESULTS/SelectedProbesets")

library(tmod)
library(dplyr)
library(tidyr)#?get
library("hgu133plus2.db")
library(stringr)
library(tibble)
library(ggplot2)
library(gplots)

#########################################################################################
# Part 1 
# - reformat data table to obtain data.frames with 2 variables: GeneID (GENE SYMBOLS from hgu133plus) and
#   value = scores. In this analysis we are only using the GeneID value.
# - eliminate probesets that are not annotated ('NA')
# - summarize genes appearing more than once in the list usign the median score.
#########################################################################################

fileNames <- list.files(path = getwd(), pattern = '.txt')
fileList <- lapply(fileNames, read.table, header = F, sep = '\t')

names(fileList) <- gsub('\\_.txt','', fileNames)

fileList <- lapply(fileList, separate, V1,into = c('geneID', 'value'), sep = ' ')
for(i in 1:length(fileList)) {
  currdata <- fileList[[i]]
  currdata <- currdata[currdata$geneID != "NA",]
  currdata$value <- as.numeric(as.character(currdata$value))
  fileList[[i]] <- currdata
}


SummarizeDuplicates <- function(df) {
  df <- df %>% group_by(geneID) %>% summarize(value=median(value))
}

summfileList <- lapply(fileList, SummarizeDuplicates)
summfileList <- lapply(summfileList, as.data.frame)

allNames <- lapply(summfileList, '[', i =, j = 1)


#########################################################################################
# Part 2 - Enrichment 
# - obtain SYMBOLS of all genes present in hgu133plus2.
# - perform an hypergeometric enrichment test using the tmod library function with a FDR
#   threshold of 0.001 (low thresholds are reccomended for the hypergeometric test in gene 
#   enrichments) and the genesets published in Li, Shuzhao, Nadine Rouphael, Sai Duraisingham, Sandra Romero-Steiner, 
#   Scott Presnell, Carl Davis, Daniel S Schmidt, et al. 2014. "Molecular Signatures of Antibody Responses
#   Derived from a Systems Biology Study of Five Human Vaccines." Nature Immunology 15 (2).
#   Nature Publishing Group: 195-204.
# - Summarize results in a bubble plot figure where the size of the bubbles is proportional to the -log10(FDR) of
#   the enrichment and the colour of the bubble is more intense when the proportion of genes in the module present in 
#   the genes in the components is high.
# - Export a table summarizing all module enrichments for each component and analytic approach, including genes that 
#   underlie the enrichment.
#########################################################################################

keys <- keys(hgu133plus2.db) 
bg <- as.character(na.omit(AnnotationDbi::select(hgu133plus2.db, keys=keys, columns = "SYMBOL")$SYMBOL))


resAll <- lapply(allNames, tmodHGtest, bg = bg, qval = 0.001)
for(i in 1:length(resAll)) {
  currdata <- resAll[[i]]
  currdata$perc <- currdata$b/currdata$B*100
  currdata$minLogP <- -log10(currdata$adj.P.Val)
  resAll[[i]] <- currdata
}

# resSumm <- lapply(resAll, '[', c('ID', 'Title', 'adj.P.Val', 'perc', 'minLogP'))
resSumm <- lapply(resAll, as.data.frame)
resSumm <- resSumm[sapply(resSumm, function(x) dim(x)[1]>0)]
for(i in 1:length(resSumm)) {
  currdata <- resSumm[[i]]
  currname <- names(resSumm)[[i]]
  currname <- substring(currname, 9)
  currdata$type <- currname
  currdata <- separate(currdata, type, into = c('analysis', '..'), sep = '_')
  resSumm[[i]] <- currdata
}
resSummTable <- do.call('rbind', resSumm)

# viz
# eliminate TBA as they get all collapsed into a single value
resSummTable2 <- filter(resSummTable, Title != 'TBA')
colnames(resSummTable2)[2] <- c('.')
resSummTable2[resSummTable2[13]=='RALL', 13] = '.'
pdf( '../annotation.pdf', width=6, height=8,paper= 'special')
ggplot(resSummTable2, aes(y = ., x = .., size=minLogP)) +        ## global aes
  geom_point(aes(colour=perc)) +    ## geom_point for circle illusion
  scale_color_gradient2(low="white", high="red" ) +
  scale_size_continuous(limits = c(1,50), range = c(0,5))+
  theme_bw() +
  # theme(plot.background = element_rect(fill="grey", colour="white"))+
  theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) +
  facet_wrap(~analysis)

dev.off()

# get genes underlying enrichments
selNames <- allNames[-5]
outputGenes <- list()
for(i in 1:length(resSumm)) {
  currdata <- resSumm[[i]]
  currgenes <- selNames[[i]]
  print(nrow(currdata))
  listGenes <- list()
  for(ii in 1:nrow(currdata)){
    currmod <- currdata$ID[ii]
    currGenes <- getGenes(currmod, genelist = currgenes )$Genes
    listGenes[[ii]] <- currGenes
  }
  outputGenes[[i]] <- listGenes
}
flatGenes <- unlist(outputGenes)

resSummTable$Genes <- flatGenes
#write.table(resSummTable, file = paste(Sys.Date(), 'WSPCAEnrichments.txt', sep = '_'), sep = '\t',
#            row.names = F)
