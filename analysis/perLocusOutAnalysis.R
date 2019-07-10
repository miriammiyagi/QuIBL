#!/usr/bin/env Rscript

suppressMessages(library("ggplot2"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("dplyr"))
suppressMessages(library("optparse"))
suppressMessages(library("reshape2"))


##### read arguments ##############

option_list = list(
  make_option(c("-i", "--inFile"), type="character", default=NULL,
              help="per locus QuIBL output file name", metavar="file path"),
  make_option(c("-o", "--outBase"), type="character", default="out",
              help="analysis output base name [default= %default]", metavar="file path")#,
  #not implemented yet
  #make_option(c("-t", "--speciesTree"), type="logical", action = "store", default=F,
  #            help="flag specifying if a species tree was given. Defaults to categorizing the most common topology as the species tree for each triplet",metavar="T/F")
);

opt_parser = OptionParser(option_list=option_list, usage = "usage: perLocusOutAnalysis.R -i inFile -o outBase [additional options]");
opt = parse_args(opt_parser);

#####read in data
tryCatch(
  {quibl <- read.csv(opt$inFile,header=T, col.names = c("tree","outgroup","branchLength","introProb"))},
  error=function(e) {
    print_help(opt_parser)
    stop("No input file supplied", call.=FALSE)
  })


outBase <- opt$outBase
taxa=unique(quibl$outgroup)
colors <- brewer.pal(8,"Paired")
quibl$outgroup <- factor(quibl$outgroup,levels=taxa)
########### make simple plots ##########

branchLengths <- ggplot() +
  geom_histogram(data=subset(quibl, outgroup==taxa[2]),aes(x=branchLength),fill=colors[2], bins=50,alpha=0.7)+
  geom_histogram(data=subset(quibl, outgroup==taxa[1]),aes(x=branchLength),fill=colors[6], bins=50,alpha=0.7)+
  geom_histogram(data=subset(quibl, outgroup==taxa[3]),aes(x=branchLength),fill=colors[4], bins=50)+
  scale_x_continuous(limits=c(0,0.07))
branchLengths
ggsave(branchLengths,file=paste0(outBase,"branchLengthsHist.pdf"),height=10,width=10)

branchLengthBox <- ggplot(data=quibl,aes(x=as.factor(outgroup),y=branchLength,col=as.factor(outgroup))) +
  geom_jitter(alpha=.1)+
  geom_boxplot(outlier.shape = NA,fill="transparent", color="black")+
  scale_color_manual(values=c(colors[6],colors[2],colors[4]), labels=taxa,name="outgroup")+
  scale_x_discrete(labels=taxa) +
  labs(x="outgroup")+
  scale_y_continuous(limits=c(0,0.07))+
  theme(legend.position = "none")
branchLengthBox
ggsave(branchLengthBox,file=paste0(outBase,"branchLengthsBox.pdf"),height=10,width=10)


introProbs <- ggplot() +
  geom_histogram(data=subset(quibl, outgroup==taxa[2]),aes(x=introProb),fill=colors[2],alpha=.7)+
  geom_histogram(data=subset(quibl, outgroup==taxa[1]),aes(x=introProb),fill=colors[6],alpha=.7)+
  geom_histogram(data=subset(quibl, outgroup==taxa[3]),aes(x=introProb),fill=colors[4])+
  labs(x="Introgression Probability",y="Number of Windows")
introProbs
ggsave(introProbs,file=paste0(outBase,"introgressionProbHist.pdf"),height=10,width=10)


BLvsIP <- ggplot(data=quibl) +
  geom_point(aes(x=branchLength,y=introProb, col=as.factor(outgroup)), alpha=.2)+
  scale_color_manual(values=c(colors[6],colors[2],colors[4]), labels=taxa)+
  scale_x_continuous(limits=c(0,0.07))+
  theme(legend.position = "none")#+
#geom_vline(xintercept =.01)
BLvsIP
ggsave(BLvsIP,file=paste0(outBase,"BranchLengthvsIntroProb.pdf"),height=10,width=10)
