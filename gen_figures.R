
rm(list=ls())

options(scipen=10)

setwd("./")

#source("./TPBayes/fastBVSR-2019v1/R-mac/bvsr.R")
source("./BicMix2/util.R")
#source("./BicMix2/gen_figures.R")


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
library("ggsci")

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

figure.path <- file.path(results.path, "figure_BicMix2")
system(paste0("mkdir -p ",figure.path))


results.path <- file.path(results.path, "results_BicMix2")
#stashDirCreate(results.path, otherrx=FALSE, grouprx=TRUE)

n.effects.list <- c(30)
#std.err.list <- c(1, 2, 3, 4, 5)
std.err.list <- 1:7

dense.list <- c("TRUE","FALSE")



gen_SFA_data3 <- function(std=2, rsd = 123, std.err=1, n.effects = 30, nfs = 20, ng = 1000, ns=200, dense=T){
    set.seed(rsd)
    
    nf <- n.effects
    nfd <- nf - nfs

    lam <- c()
    lams <- c()
    lamd <- c()

    lams <- matrix(0,nrow=ng,ncol=nfs)

    ########## simulate lam
    block <- 5
    for(i in 1:nfs){
        #start <- sample(1:(ng-block),1)
        n.nonzero <- block + sample(1:5,1)
        index <- sample(1:ng,n.nonzero,replace=T)
        lams[index,i] = rnorm(n.nonzero,2,std)
        lams[sample(index,n.nonzero/2),i] = rnorm(n.nonzero/2,-2,std)

    }
    
    if(dense){
        l <- ng*nfd
        lamd.tmp <- rnorm(l,2,std)
        lamd.tmp[sample(1:l,l/2)] <- rnorm(l/2,-2,std)
        lamd <- matrix(lamd.tmp,nrow=ng,ncol=nfd)
        lam <- cbind(lamd,lams)
        #lam <- lam[,sample(1:nf,nf,replace=F)]
    }else{
        lam <- lams
        #lam <- lam[,sample(1:nfs,nfs,replace=F)]
    }

    ############### simulate ex
    if(!dense){
        nf <- nfs
    }
    
    ex <- matrix(rnorm(ns*nf,0,std),nrow=nf,ncol=ns)
   
    ###################### simulate error
    err <- matrix(rnorm(ng*ns,0,std.err),nrow=ng,ncol=ns)

    y <- lam %*% ex + err
    
    if(dense){
        return(list(y=y,lams=lams,lamd=lamd,lam=lam,ex=ex, err=err))
    }else{
        return(list(y=y,lams=lams,lam=lam,ex=ex, err=err))
    }
}


mine_heatmap <- function(x,title=NA){
    library(reshape2)
    library("ggsci")
    m <- melt(x)
    p <- ggplot(m,aes(x=Var1,y=Var2)) + geom_tile(aes(fill=value),color="black") + coord_flip() +theme_minimal() + theme(legend.position="none",axis.title = element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),plot.title = element_text(face = "bold", hjust = 0.5,vjust=-2),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(title) + scale_fill_gradient2(low="blue",mid="white",high="red") 
    #scale_fill_gradient2(low="green",mid="white",high="red") 
    #panel.border = element_rect(colour = "black", fill=NA, size=1)
}

sfa_scheme <- function(){
    #data=gen_SFA_data2(std=2, rsd = 123, std.err=1, n.effects = 5, nfs=3, ng=10,ns=5,dense=TRUE)
    data=gen_SFA_data3(std=0.5, rsd = 123, std.err=1, n.effects = 10, nfs=7, ng=20,ns=20,dense=TRUE)

    plam <- mine_heatmap(data$lam,expression(Lambda["P x K"]))
    pex <- mine_heatmap(data$ex,expression(X["K x N"]))
    py <- mine_heatmap(data$y,expression(Y["P x N"]))
    perr <- mine_heatmap(data$err,expression(bold(epsilon)["P x N"]))
    return(list(plam = plam, pex = pex, py = py, perr = perr))
}




p <- sfa_scheme()

png(file.path(figure.path,"Lam.png"),width=200,height=800,res=300)
print(p$plam)
dev.off()


png(file.path(figure.path,"EX.png"),width=800,height=400,res=300)
print(p$pex)
dev.off()

png(file.path(figure.path,"Y.png"),width=800,height=800,res=300)
print(p$py)
dev.off()

png(file.path(figure.path,"Err.png"),width=800,height=800,res=300)
print(p$perr)
dev.off()


height=1
width=2
height2 = 0.5
pdf(file.path(figure.path,"Lam.pdf"),width=height,height=width)
print(p$plam)
dev.off()


pdf(file.path(figure.path,"EX.pdf"),width=width,height=height)
print(p$pex)
dev.off()

pdf(file.path(figure.path,"Y.pdf"),width=width,height=width)
print(p$py)
dev.off()

pdf(file.path(figure.path,"err.pdf"),width=width,height=width)
print(p$perr)
dev.off()





