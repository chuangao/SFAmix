

rm(list=ls())

options(scipen=10)

setwd("./")

#source("./TPBayes/fastBVSR-2019v1/R-mac/bvsr.R")
source("./BicMix2/util.R")

set.seed(123)
library(ggplot2)
library(TPBayes)
library(parallel)
library(glmnet)
library(gridExtra)
library(SSLASSO)
library(dplyr)
library(ggpubr)

library(BicMix)
library(BicMix2)
library(PMA)
numCores <- 6

#library(knitr)
#require(markdown)

#setwd("./code/R")

## file directories
results.path <- "."

## testing directory creations
table.path <- file.path(results.path, "table_BicMix2")
#stashDirCreate(table.path, otherrx=FALSE, grouprx=TRUE)

plot.path <- file.path(results.path, "plot_BicMix2")
#stashDirCreate(plot.path, otherrx=FALSE, grouprx=TRUE)
data.path <- file.path(results.path, "data_BicMix2")

results.path <- file.path(results.path, "results_BicMix2")
#stashDirCreate(results.path, otherrx=FALSE, grouprx=TRUE)

n.effects.list <- c(30)
#std.err.list <- c(1, 2, 3, 4, 5)
std.err.list <- 1:7

dense.list <- c("TRUE","FALSE")

method.list <- c("BicMix2","BicMix","IFA","BFRM","SPCA","KSVD")
method.list <- c("BicMix2")

sed.list <- 1:5
std.effect.list <- 2
b.list <- c(10,100,1000, 10000) * 10000

#b.list <- c(10000)*10000
sn <- 500
dy <- 200

ng=sn
ns=dy

nfs = 40
nf = 50

a=0.1

min_fac <- 20
step <- 2
#b=1000000
#c=0.1; d=1000; g=0.1; h=100000


param.config <- expand.grid(n.effects.list, std.err.list, method.list, sed.list, std.effect.list, b.list,dense.list)
names(param.config) <- c("n.effects", "std.err", "method", "seed", "std.effect","b","dense")

#gen_bulk_data(param.config,output.dir=data.path)

#ksvd <- read.table("results_BicMix2/KSVD/n.effects30_std.err7_seed10_std.effect2_b1000000.coef")
#plot(as.numeric(ksvd[1,]))

source("./BicMix2/KSVD/ksvd.r")
source("./BicMix2/IFA/ifa.r")
source("./BicMix2/BFRM/bfrm.r")
source("./BicMix2/util.R")

method <- "BicMix2"

itr <- 1001
#res <- run_sim(param.config[param.config$method == method & param.config$dense == TRUE,][1,], nfs=nfs, nf=nf, ng= ng, ns=ns, itr = itr)
res <- run_sim(param.config, itr = itr, nfs=nfs, nf=nf, ng=sn, ns=dy,min_fac=min_fac,step=step)
#results <- res[[1]]
#count.prob <- apply(results$z,2,function(x){return(sum(x > 0.5)/sn)})

for(i in 1:length(res)){
    names(res[[i]]$score.dense)[names(res[[i]]$score.dense) == "score.dense.1"] <- "score.dense.precis" 
}

res2 <- extract_res(res)
#res2.bak <- res2
library("ggsci")

score.sparse.all <- res2$score.sparse
score.dense.all <- res2$score.dense



#write.csv(score.sparse.all,file.path(table.path,"comparison_only_sparse_stdErr1_5_stdEff2_p500_n200_a0.5_b1000000.csv"),row.names=F)
#write.csv(score.dense.all,file.path(table.path,"comparison_sparse_with_dense_stdErr1_5_stdEff2_p500_n200_a0.5_b1000000.csv"),row.names=F)
#score.sparse.all <- read.csv(file.path(table.path,"comparison_only_sparse_stdErr1_5_stdEff2_p500_n200_a0.5_b10000.csv"))
#score.dense.all <- read.csv(file.path(table.path,"comparison_sparse_with_dense_stdErr1_5_stdEff2_p500_n200_a0.5_b10000.csv"))



score.sparse.all$std.err <- factor(score.sparse.all$std.err,levels=unique(score.sparse.all$std.err))
score.dense.all$std.err <- factor(score.dense.all$std.err,levels=unique(score.dense.all$std.err))

score.sparse.all$method <- factor(score.sparse.all$method,levels=unique(score.sparse.all$method))
score.dense.all$method <- factor(score.dense.all$method,levels=unique(score.dense.all$method))

names(score.sparse.all) <- c("n.effects", "Std.err", "seed",  "Method", "Std.effect","b", "dense", "Score","Score.precis","nfs","nfd")
names(score.dense.all) <- c("n.effects", "Std.err", "seed",  "Method", "Std.effect","b", "dense", "Score","Score.precis","nfs","nfd")

score.dense.all <- score.dense.all[score.dense.all$dense == TRUE,]

score.sparse.all$Data <- ifelse(score.sparse.all$dense,"Sparse-sim1","Sparse-sim2")
score.dense.all$Data <- ifelse(score.dense.all$dense,"Dense-sim1")

score.all <- rbind(score.sparse.all,score.dense.all)
score.all$Data <- factor(score.all$Data,levels=unique(score.all$Data))

#file.csv <- file.path(table.path,paste0("comparison_self_stdErr1_7_stdEff2_p",sn,"_n",dy,"_nf",nf,"_nfs",nfs,"_nf_start100_a0.1_bstart1_5_10_50_100_1000E4_nflt1_bhalf_bfGrt1000.csv"))
file.csv <- file.path(table.path,paste0("comparison_self_stdErr1_7_stdEff2_p",sn,"_n",dy,"_nf",nf,"_nfs",nfs,"_nf_start100_a",0.1,"_bstart10_100_1000_10000E4_min_fac",min_fac,"_bstep",step,"_zGrt0.5.csv"))

#write.csv(score.all,file.csv,row.names=F)

score.all <- read.csv(file.csv,as.is=T) 


#score.all$b <- factor(score.all$b,levels=unique(score.all$b))
score.all$Std.err <- factor(score.all$Std.err,levels=unique(score.all$Std.err))
#score.all$Std.err <- paste0("Sd.err=",score.all$Std.err)
score.all$Method[score.all$Method == "BicMix2"] <- "BARF"
#score.all$Method <- factor(score.all$Method,levels=c("BARF","BicMix","BFRM","IFA","KSVD","SPCA"))
score.all$Data <- factor(score.all$Data,levels=c("Sparse-sim1","Sparse-sim2","Dense-sim1"))
score.all$b <- factor(paste0("b=",score.all$b),levels=unique(paste0("b=",score.all$b)))
#p <- ggplot(score.all, aes(x=Std.err,y=Score, color=Method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="top",legend.title=element_blank()) + scale_color_jco() + facet_grid(Data ~ . , scale="free_y")

p <- ggplot(score.all, aes(x=Std.err,y=Score)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="top",legend.title=element_blank()) + scale_color_jco() + facet_grid(Data ~ b , scale="free_y") + xlab("Sd.err") + ylim(0,1)

file.pdf <- file.path(plot.path,paste0("comparison_self_stdErr1_7_stdEff2_p",sn,"_n",dy,"_nf",nf,"_nfs",nfs,"_nf_start100_a",0.1,"_bstart10_100_1000_10000E4_min_fac",min_fac,"_bstep",step,"_zGrt0.5.pdf"))
pdf(file.pdf,width=6,height=4)
print(p)
dev.off()


p <- ggplot(score.all, aes(x=Std.err,y=Score.precis)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="top",legend.title=element_blank()) + scale_color_jco() + facet_grid(Data ~ b , scale="free_y") + xlab("Sd.err") +ylim(0,1)
 
file.pdf <- file.path(plot.path,paste0("comparison_self_stdErr1_7_stdEff2_p",sn,"_n",dy,"_nf",nf,"_nfs",nfs,"_nf_start100_a",0.1,"_bstart10_100_1000_10000E4_min_fac",min_fac,"_bstep",step,"_zGrt0.5_precis.pdf"))
pdf(file.pdf,width=6,height=4)
print(p)
dev.off()


score.all$nf <- score.all$nfd + score.all$nfs
p.nf <- ggplot(score.all, aes(x=Std.err,y=nf)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="top",legend.title=element_blank()) + scale_color_jco() + facet_grid(Data ~ b , scale="free_y") + xlab("Sd.err") + ylab("N Factors")


file.pdf <- file.path(plot.path,paste0("nf_comparison_self_stdErr1_7_stdEff2_p",sn,"_n",dy,"_nf",nf,"_nfs",nfs,"_nf_start100_a",0.1,"_bstart10_100_1000_10000E4_min_fac",min_fac,"_bstep",step,"_zGrt0.5.pdf"))
pdf(file.pdf,width=6,height=4)
print(p.nf)
dev.off()
























score.sparse.all$std.err <- factor(score.sparse.all$std.err,levels=unique(score.sparse.all$std.err))
score.dense.all$std.err <- factor(score.dense.all$std.err,levels=unique(score.dense.all$std.err))

score.sparse.all$method <- factor(score.sparse.all$method,levels=unique(score.sparse.all$method))
score.dense.all$method <- factor(score.dense.all$method,levels=unique(score.dense.all$method))

score.sparse.all

############# plot up results with different b values
p.sparse1 <- ggplot(score.sparse.all[score.sparse.all$dense == TRUE,], aes(x=std.err,y=score.sparse, color=method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="top",legend.title=element_blank()) + scale_color_jco() + ylim(0,1) + facet_grid(. ~ b)
#ggthemes::geom_tufteboxplot()
p.sparse2 <- ggplot(score.sparse.all[score.sparse.all$dense == FALSE,], aes(x=std.err,y=score.sparse, color=method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="none",legend.title=element_blank()) + scale_color_jco() + ylim(0,1) + facet_grid(. ~ b)

lg <- as_ggplot(get_legend(p.sparse1)) + theme(plot.margin = unit(c(0,0,0,0), "cm"))
   
p.sparse1 <- p.sparse1 + theme(legend.position="none")

p.dense <- ggplot(score.dense.all[score.sparse.all$dense == TRUE,], aes(x=std.err,y=score.dense, color=method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="none") + scale_color_jco() + facet_grid(. ~ b)


#lay <- matrix(c(rep(c(1,2),6),3,3),nrow=2)
#lay <- matrix(rep(1:4,c(1,4,4,4)),ncol=1)
#lay <- matrix(c(c(1,2,2,2,2),c(1,3,3,3,3),c(1,4,4,4,4)),ncol=3)
lay <- matrix(c(
    rep(1:3,6),
    rep(4,3)
),nrow=3)


file.png <- file.path(plot.path,paste0("comparison_self_stdErr1_7_stdEff2_p",sn,"_n",dy,"_nf",nf,"_nfs",nfs,"_a0.1_bstart1_5_10_50_100E4_nflt1_end_bfGrt1000.png"))
png(file.png,width=2400,height=1600, res = 300)
grid.arrange(p.sparse1,p.sparse2,p.dense,ncol=1) + theme(plot.margin = unit(c(0,0,0,0), "cm"))
dev.off()


file.pdf <- file.path(plot.path,paste0("comparison_self_stdErr1_7_stdEff2_p",sn,"_n",dy,"_nf",nf,"_nfs",nfs,"_a0.1_bstart1_5_10_50_100E4_nflt1_end_bfGrt1000.pdf"))
pdf(file.pdf,width=8,height=6)
grid.arrange(p.sparse1,p.sparse2,p.dense,ncol=1) + theme(plot.margin = unit(c(0,0,0,0), "cm"))
dev.off()


############# plot up results with different b values

score.sparse.all$b <- factor(score.sparse.all$b,levels=unique(score.sparse.all$b))
score.dense.all$b <- factor(score.dense.all$b,levels=unique(score.dense.all$b))

p.sparse1 <- ggplot(score.sparse.all[score.sparse.all$dense == TRUE,], aes(x=std.err,y=score.sparse, color=b)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="top",legend.title=element_blank()) + scale_color_jco() + ylim(0,1) 
#ggthemes::geom_tufteboxplot()
p.sparse2 <- ggplot(score.sparse.all[score.sparse.all$dense == FALSE,], aes(x=std.err,y=score.sparse, color=b)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="none",legend.title=element_blank()) + scale_color_jco() + ylim(0,1)

lg <- as_ggplot(get_legend(p.sparse1)) + theme(plot.margin = unit(c(0,0,0,0), "cm"))
   
p.sparse1 <- p.sparse1 + theme(legend.position="none")

p.dense <- ggplot(score.dense.all[score.sparse.all$dense == TRUE,], aes(x=std.err,y=score.dense, color=b)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="none") + scale_color_jco()


#lay <- matrix(c(rep(c(1,2),6),3,3),nrow=2)
#lay <- matrix(rep(1:4,c(1,4,4,4)),ncol=1)
lay <- matrix(c(c(1,2,2,2,2),c(1,3,3,3,3),c(1,4,4,4,4)),ncol=3)

grid.arrange(lg, p.sparse1,p.sparse2,p.dense, layout_matrix = lay) + theme(plot.margin = unit(c(0,0,0,0), "cm"))
