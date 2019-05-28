#!/usr/bin/env Rscript


library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
library("optparse")

##### read arguments ##############

option_list = list(
  make_option(c("-i", "--inFile"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--outBase"), type="character", default="out",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-s", "--sigLevel"), type="integer", default=-10,
              help="delta BIC value for significance [default= %default]", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



quibl <- read.csv("~/Desktop/testOut.csv",header=T,
                  colClasses = c(rep("factor",2),rep("numeric",9)))

#####read in data
quibl <- read.csv(opt$inFile,header=T,
                  colClasses = c(rep("factor",2),rep("numeric",9)))

outBase <- opt$outBase
sigLevel <- opt$sigLevel

## add columns for delta BIC, total introgression proportion, significance, most common tree in triplet (if species tree not specified)
quibl$BICdiff=quibl$BIC1-quibl$BIC2

#get the total loci analyzed by summing all values of numtrees for a single triplet set
numLoci <- sum(subset(quibl,triplet==quibl$triplet[1])$count)
quibl$totalIntroProp <- quibl$mixprop2*(quibl$count/numLoci)
quibl$isSig <- quibl$BICdiff < sigLevel
quibl$mostCommon <- F
for (trip in unique(as.character(quibl$triplet))){
  maxVal=max(subset(quibl,triplet==trip)$count)
  quibl[which(quibl$triplet==trip & quibl$count==maxVal),]$mostCommon <-  T
}


colors=brewer.pal(8,"Paired")

##### total introgression proportions #####
#here, we're calculating the total number of discordant triplets that arose via introgression. Assuming that concordant triplets are most common.
discordSet=subset(quibl, mostCommon==F)
allDiscordTrees=sum(discordSet$count)
pctDiscordTrees=allDiscordTrees/sum(quibl$count)
discordIntro <- sum(discordSet$count*discordSet$mixprop2*as.numeric(discordSet$isSig))
propDiscordIntro <- discordIntro/allDiscordTrees
propOverallIntro <- discordIntro/sum(quibl$count)
print(paste0("Number of discordant topologies with significant evidence for introgression: ",length(which(discordSet$isSig==T))))
print(paste0("Proportion of discordant trees arising via introgression: ",propDiscordIntro))
print(paste0("Proportion of all trees arising via introgression: ",propOverallIntro))

#### overview graphs
CvalHist <- ggplot(data=quibl,aes(x=C2))+
  geom_histogram(bins=20)+
  labs(y="Number of topologies",x="Non-ILS C")
#geom_vline(xintercept = 0.5)
CvalHist
ggsave(CvalHist,file=paste0(outBase,"_CVal_hist.pdf"),height=10,width=10)

introPropHist <- ggplot(data=quibl,aes(x=mixprop2))+
  geom_histogram(bins=20)+
  labs(y="Number of topologies",x="Topology Non-ILS proportion")
introPropHist
ggsave(introPropHist,file=paste0(outBase,"_topIntroProp_hist.pdf"),height=10,width=10)

CvalVsIntroProp <- ggplot(data=quibl)+
  geom_point(aes(x=C2,y=mixprop2,size=count/numLoci, col=mostCommon, pch=isSig), alpha=.7)+
  labs(y="Topology non-ILS proportion",x="Non-ILS C")+
  scale_size_continuous(name = "Proportion of topologies")+
  scale_color_manual(values=c(colors[8],colors[2]),labels=c("minor topology","major topology"),name="")+
  scale_shape_manual(values=c(18,19),labels=c("nonsignificant","significant"),name="")#+
  #theme(legend.position="none")
CvalVsIntroProp
ggsave(CvalVsIntroProp,file=paste0(out,"CValVsTopIntroProp.pdf"),height=10,width=10)

CvalVsTotalIntroProp <- ggplot(data=quibl)+
  geom_point(aes(x=C2,y=totalIntroProp,size=count/numLoci, col=mostCommon, pch=isSig), alpha=.7)+
  labs(y="Total Non-ILS proportion",x="Non-ILS C")+
  scale_size_continuous(name = "Proportion of topologies")+
  scale_color_manual(values=c(colors[8],colors[2]),labels=c("minor topology","major topology"),name="")+
  scale_shape_manual(values=c(18,19),labels=c("nonsignificant","significant"),name="")#+
  #theme(legend.position="none")
CvalVsTotalIntroProp
ggsave(CvalVsTotalIntroProp,file=paste0(out,"CValVsTotalIntroProp.pdf"),height=10,width=10)
