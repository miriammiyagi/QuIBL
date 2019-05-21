### #!/usr/bin/env Rscript


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
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);





#####read in data
#currently, output doesn't have numTrees. when this is fixed, sub in.
#quibl <- read.csv(opt$inFile,header=F,col.names=c("triplet","outgroup","C1","C2","prop1","prop2","numTrees","BIC1","BIC2"))
quibl <- read.csv(opt$inFile,header=F,col.names=c("triplet","outgroup","C1","C2","prop1","prop2","BIC1","BIC2"))

## add columns for delta BIC, total introgression proportion, significance, most common tree in triplet (if species tree not specified)
quibl$BICdiff=quibl$BIC1-quibl$BIC2

#get the total loci analyzed by summing all values of numtrees for a single triplet set
numLoci <- sum(subset(quibl,triplet==quibl$triplet[1])$numTrees)
#quibl$totalIntroProp <- quibl$prop2*(quibl$numTrees/numLoci)
quibl$isSig <- quibl$BICdiff < (-10)
quibl$mostCommon <- F
quibl$mostCommon[match(unique(quibl$triplet), quibl$triplet)] <- T


colors=brewer.pal(8,"Paired")

##### total introgression proportions #####
#here, we're calculating the total number of discordant triplets that arose via introgression. Assuming that concordant triplets are most common.
discordSet=subset(quibl, mostCommon==F)
allDiscordTrees=sum(discordSet$numTrees)
pctDiscordTrees=allDiscordTrees/sum(quibl$numTrees)
#discordIntro <- sum(discordSet$numTrees*discordSet$prop2*as.numeric(discordSet$isSig))
discordIntro <- sum(discordSet$numTrees*discordSet$prop2)
propDiscordIntro <- discordIntro/allDiscordTrees
propOverallIntro <- discordIntro/sum(quibl$numTrees)


#### overview graphs
CvalHist <- ggplot(data=quibl,aes(x=C2))+
  geom_histogram(bins=20)+
  labs(y="Number of topologies",x="introgression C")
  #geom_vline(xintercept = 0.5)
CvalHist
ggsave(CvalHist,file=paste0(outBase,"_CVal_hist.pdf"),height=10,width=10)

introPropHist <- ggplot(data=quibl,aes(x=prop2))+
  geom_histogram(bins=20)+
  labs(y="Number of topologies",x="Topology introgression proportion")
introPropHist
ggsave(introPropHist,file=paste0(outBase,"_topIntroProp_hist.pdf"),height=10,width=10)

CvalVsIntroProp <- ggplot(data=quibl)+
  geom_point(aes(x=C2,y=prop2,size=numTrees/numLoci, col=mostCommon, pch=isSig), alpha=.7)+
  labs(y="Topology introgression proportion",x="Introgression C")+
  scale_size_continuous(name = "Proportion of topologies")+
  scale_color_manual(values=c(colors[8],colors[2]),labels=c("minority","majority"),name="")+
  scale_shape_manual(values=c(18,19),labels=c("nonsignificant","significant"),name="")+
  theme(legend.position="none")
CvalVsIntroProp
ggsave(CvalVsIntroProp,file=paste0(out,"CValVsTopIntroProp.pdf"),height=10,width=10)

CvalVsTotalIntroProp <- ggplot(data=combined)+
  geom_point(aes(x=C2,y=totalIntroProp,size=numTrees/numLoci, col=mostCommon, pch=isSig), alpha=.7)+
  labs(y="Total introgression proportion",x="Introgression C")+
  scale_size_continuous(name = "Proportion of topologies")+
  scale_color_manual(values=c(colors[8],colors[2]),labels=c("minority","majority"),name="")+
  scale_shape_manual(values=c(18,19),labels=c("nonsignificant","significant"),name="")+
  theme(legend.position="none")
CvalVsTotalIntroProp
ggsave(CvalVsTotalIntroProp,file=paste0(out,"CValVsTotalIntroProp.pdf"),height=10,width=10)
