


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
std.err.list <- 1:5

dense.list <- c("TRUE","FALSE")

method.list <- c("BicMix2","BicMix")
sed.list <- 1:5
std.effect.list <- 2
b.list <- c(1000000)
sn <- 500
ng <- 1000 

a=0.1
#b=1000000
#c=0.1; d=1000; g=0.1; h=100000


param.config <- expand.grid(n.effects.list, std.err.list, method.list, sed.list, std.effect.list, b.list,dense.list,sn.list)
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
res <- run_sim(param.config[param.config$method == method & param.config$dense == FALSE,][1,], itr = itr)
#res <- run_sim(param.config, itr = itr)
#results <- res[[1]]
#count.prob <- apply(results$z,2,function(x){return(sum(x > 0.5)/sn)})

res2 <- extract_res(res)
#res2.bak <- res2
library("ggsci")

score.sparse.all <- res2$score.sparse
score.dense.all <- res2$score.dense

#write.csv(score.sparse.all,file.path(table.path,"comparison_only_sparse_stdErr1_5_stdEff2_p500_n200_a0.5_b10000.csv"),row.names=F)
#write.csv(score.dense.all,file.path(table.path,"comparison_sparse_with_dense_stdErr1_5_stdEff2_p500_n200_a0.5_b10000.csv"),row.names=F)
#score.sparse.all <- read.csv(file.path(table.path,"comparison_only_sparse_stdErr1_5_stdEff2_p500_n200_a0.5_b10000.csv"))
#score.dense.all <- read.csv(file.path(table.path,"comparison_sparse_with_dense_stdErr1_5_stdEff2_p500_n200_a0.5_b10000.csv"))


# score.sparse.all <- rbind(score.sparse1, score.sparse)
# score.dense.all <- rbind(score.dense1, score.dense)

# write.csv(score.sparse.all,file.path(table.path,paste0(method,"_score_sparse.csv")),row.names=F)
# #write.csv(score.sparse.all,file.path(table.path,paste0(method,"_spca_ksvd_score_sparse.csv")),row.names=F)

# write.csv(score.dense.all,file.path(table.path,paste0(method,"_score_dense.csv")),row.names=F)

# score.sparse.all <- do.call(rbind,lapply(method.list,function(x){
#     return(read.csv(file.path(table.path,paste0(x,"_score_sparse.csv"))))
# }))

# score.dense.all <- do.call(rbind,lapply(method.list,function(x){
#     return(read.csv(file.path(table.path,paste0(x,"_score_dense.csv"))))
# }))

score.sparse.all$std.err <- factor(score.sparse.all$std.err,levels=unique(score.sparse.all$std.err))
score.dense.all$std.err <- factor(score.dense.all$std.err,levels=unique(score.dense.all$std.err))

score.sparse.all$method <- factor(score.sparse.all$method,levels=unique(score.sparse.all$method))
score.dense.all$method <- factor(score.dense.all$method,levels=unique(score.dense.all$method))

score.all 
########### plot up the score for data with dense components
p.sparse1 <- ggplot(score.sparse.all[score.sparse.all$dense == TRUE,], aes(x=std.err,y=score.sparse, color=method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="top",legend.title=element_blank()) + scale_color_jco() + ylim(0,1)
#ggthemes::geom_tufteboxplot()
p.sparse2 <- ggplot(score.sparse.all[score.sparse.all$dense == FALSE,], aes(x=std.err,y=score.sparse, color=method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="none",legend.title=element_blank()) + scale_color_jco() + ylim(0,1)

lg <- as_ggplot(get_legend(p.sparse1)) + theme(plot.margin = unit(c(0,0,0,0), "cm"))
   
p.sparse1 <- p.sparse1 + theme(legend.position="none")

p.dense <- ggplot(score.dense.all[score.sparse.all$dense == TRUE,], aes(x=std.err,y=score.dense, color=method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="none") + scale_color_jco() + ylim(0,1)


#lay <- matrix(c(rep(c(1,2),6),3,3),nrow=2)
#lay <- matrix(rep(1:4,c(1,4,4,4)),ncol=1)
lay <- matrix(c(c(1,2,2,2,2),c(1,3,3,3,3),c(1,4,4,4,4)),ncol=3)


file.png <- file.path(plot.path,paste0("comparison_stdErr1_5_stdEff2_p500_n200_a0.5_b10000.png"))
png(file.png,width=2400,height=1600, res = 300)
grid.arrange(lg, p.sparse1,p.sparse2,p.dense,layout_matrix = lay) + theme(plot.margin = unit(c(0,0,0,0), "cm"))
dev.off()


########### plot up the score for data with only components

p.sparse2 <- ggplot(score.sparse.all[score.sparse.all$dense == FALSE,], aes(x=std.err,y=score.sparse, colour=method)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="right") + scale_color_jco() + ylim(0,1)

file.png <- file.path(plot.path,paste0("BicMix2_SPCA_KSVD_1atom_score_only_sparse.png"))
png(file.png,width=1600,height=1200, res = 300)
print(p.sparse2)
dev.off()




        i=1
        n.effects <- param.config[i,"n.effects"]
        std.err <- param.config[i,"std.err"]
        i.seed <- param.config[i,"seed"]
        method <- param.config[i,"method"]
        std.effect <- param.config[i,"std.effect"]
        b <- param.config[i,"b"]
        

        a=0.1; b=1000000
        std.err = 1
        #c=0.1; d=10000; g=0.1; h=10000


        set.seed(123 * i.seed)
        
        data=gen_SFA_data(std=2, rsd = 123, std.err=std.err, n.effects = 30, ng=1000,ns=200)
        
        #results.dir <- results.path
        #data.dir <- "data"
        
        #y <- t(apply(data$y,1,function(x){return(qqnorm(x,plot=F)$x)}))
        y <- data$y
        results <- BicMixR2(y, nf = 100, a=a,b=b,c=c,d=d,g=g,h=h,rsd=123,itr=2001, tol = 0.001, out_dir=results.path, out_itr = 100)

        #apply(cor(data$lams),2,function(x){return(max(abs(x[x!= 1])))})

        count.prob <- apply(results$z,2,function(x){return(sum(x > 0.5)/sn)})
        lams <- results$lam[,count.prob < 0.5]
        lamd <- results$lam[,count.prob >= 0.5]
        #corr <- cor(results$lam[,count.prob < 0.5 & count.prob > 0.01],data$lams)
        corr.sparse <- cor(lams,data$lams)

        corr.dense <- cor(lamd,data$lamd)


        image(corr.sparse)
        #image(corr.dense)

        cal_score_sparse(corr.sparse)
        cal_score_dense(lamd,data$lamd)


        #var.comp <- apply(results$lam,2,function(x){return(sd(x))})

        library(reshape2)
        library(ggplot2)

        lam.melt <- melt(results$lam)
        names(lam.melt) <- c("Index","Factor","Value")
        p.lam <- ggplot(lam.melt,aes(x=Index,y=Value)) + geom_point() + facet_grid(Factor ~ .)

        z.melt <- melt(results$z)
        names(z.melt) <- c("Index","Factor","Value")
        p.z <- ggplot(z.melt,aes(x=Index,y=Value)) + geom_point() + facet_grid(Factor ~ .)


        png(file.path(plot.path,"factor_scatter_plot_recovered.png"),width=1200,height=16000)
        grid.arrange(p.lam,p.z,ncol=2)
        dev.off()
















beta.all <- c()
pred.all <- c()
recover.all <- c()

itr <- 501

res <- run_sim(param.config[1:12,], itr = itr, output.dir = results.path)
res2 <- extract_res(res)
beta.all <- rbind(beta.all, res2$betas)
pred.all <- rbind(pred.all, res2$pred)
recover.all <- rbind(recover.all, res2$recover)

# file <- file.path(table.path,paste0("betas_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_to_others_itr",itr,".csv"))
# write.csv(beta.all,file,row.names=F)
# file <- file.path(table.path,paste0("pred_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_to_others_itr",itr,".csv"))
# write.csv(pred.all,file,row.names=F)
# file <- file.path(table.path,paste0("recover_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_to_others_itr",itr,".csv"))
# write.csv(recover.all,file,row.names=F)


beta.all$prob.true <- ifelse(beta.all$beta == 0, 0, 1)

beta.all.tpb.cut <- beta.all[beta.all$method == "TPB",]
beta.all.tpb.cut$prob <- ifelse(beta.all.tpb.cut$prob > 0.5, 1,0)
beta.all.tpb.cut$method <- rep("TPB_cut",nrow(beta.all.tpb.cut))

beta.all.bvsr.cut <- beta.all[beta.all$method == "BVSR",]
beta.all.bvsr.cut$prob <- ifelse(beta.all.bvsr.cut$prob > 0.5, 1,0)
beta.all.bvsr.cut$method <- rep("BVSR_cut",nrow(beta.all.bvsr.cut))

beta.all <- rbind(beta.all, beta.all.tpb.cut)
beta.all <- rbind(beta.all, beta.all.bvsr.cut)

#beta.all %>% filter(method=="TPB") %>% group_by(n.effects,std.err,seed,shape) %>% summarise(sum(prob > 0.001))

library("plotROC")
library("pROC")
#library(ROCit)
#library(ROCR)
library("ggsci")
#library(viridis)


beta.roc <- beta.all %>% group_by(n.effects,std.err,seed,method,shape) %>% do(cal_roc(.))
beta.roc2 <- beta.roc %>% group_by(n.effects,std.err,seed,method,shape) %>% do(cal_quant(df=.))

#with(beta.all[1:10000,],sum(prob > 0.5 & prob.true == 1))/20

#beta.roc %>% group_by(n.effects,std.err,seed,method,shape) %>% summarise(n())

beta.roc2 <- beta.roc2 %>% group_by(n.effects,std.err,method,shape,fp) %>% summarise(tpm=mean(tp),tpmin=min(tp),tpmax=max(tp)) 
beta.roc2$std.err <- paste0("Sd.err=",beta.roc2$std.err)
beta.roc2$std.err <- factor(beta.roc2$std.err,levels=unique(beta.roc2$std.err))

beta.roc2$n.effects <- paste0("N.effects=",beta.roc2$n.effects)
beta.roc2$n.effects <- factor(beta.roc2$n.effects,levels=unique(beta.roc2$n.effects))

beta.roc2$Methods <-  beta.roc2$method
#beta.roc2$Binary <- factor(ifelse(beta.roc2$Methods == "TPB_cut" | beta.roc2$Methods == "BVSR_cut" | beta.roc2$Methods == "SSLASSO", "Yes", "No"),levels=c("No","Yes"))
for(i in 1:length(shape.list)){
    shape <- shape.list[i]

    p <- ggplot(beta.roc2[beta.roc2$shape == shape,], aes(x=fp,y=tpm)) + geom_line(aes(colour=Methods,linetype=Methods)) + geom_ribbon(aes(ymin=tpmin,ymax=tpmax,fill=Methods),alpha=0.2,colour = NA) + facet_grid(n.effects ~ std.err) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="top") + xlab("False Positive Rate") + ylab("True Positive Rate") + scale_fill_jco() + scale_color_jco()
   
    file.png <- file.path(plot.path,paste0("roc_shape_",shape,"_nf",nf,"_comparing_method_",method,"_a",a,"_d_",d,"_to_others_itr",itr,".png"))
    png(file.png,width=2400,height=2000, res = 300)
    print(p)
    dev.off()
}




























































####################################################
if(1==0){
for(i in 1:length(std.err.list)){
    for(j in 1:length(shape.list)){
        std.err <- std.err.list[i]

        ## plot betas
        p <- ggplot(beta.all[beta.all$std.err ==std.err,], aes(lam,beta)) + geom_point(shape=1) + facet_grid(method ~ n.effects) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        file.png <- paste0("plots/betas_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_sd_err",std.err,"_to_others.png")
        png(file.png,width=1600,height=1200, res = 300)
        print(p)
        dev.off()

        ## plot predictions
        p <- ggplot(pred.all[pred.all$std.err ==std.err,], aes(y.true,y.pred)) + geom_point(shape=1) + facet_grid(method ~ n.effects) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        file.png <- paste0("plots/pred_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_sd_err",std.err,"_to_others.png")

        png(file.png,width=1600,height=1200, res = 300)
        print(p)
        dev.off()
    }
}

sum(beta.all[beta.all$method=="SSLASSO",][1:1000,]$prob)
sum(beta.all[beta.all$method=="LASSO",][1:1000,]$prob)

for(i in 1:length(shape.list)){
    shape <- shape.list[i]
    p <- ggplot(recover.all[recover.all$shape == shape,], aes(y=Recover,x=method)) + geom_boxplot(width=0.5, outlier.shape = 20, outlier.size=1) + ylim(0,1) +
     facet_grid(std.err ~ n.effects) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

}

beta.i <- beta.all[1:1000,]
rocobj <- roc(beta.i$prob.true, beta.i$prob)

file.csv <- paste0("tables/betas_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_to_others.csv")
beta <- read.csv(file.csv)




















recover <- do.call(rbind,lapply(((1:nrow(param.config))),function(i){
    cat(i,"\n")

    n.effects <- param.config[i,"n.effects"]
    std.err <- param.config[i,"std.err"]
    i.seed <- param.config[i,"seed"]
    
    beta.i <- beta[beta$n.effects==n.effects & beta$std.err==std.err & beta$seed==i.seed,]
    
    ntotal <- sum(beta.i[beta.i$method == "TPBayes",]$beta !=0 )

    #rtpb <- sum(beta.i[beta.i$method == "TPBayes",]$prob > 0.99)/ntotal
    obj <- kmeans(abs(beta.i$lam[beta.i$method == "TPBayes"]),2)
    rtpb <- sum(obj$cluster == which.max(obj$centers))/ntotal
    
    obj <- kmeans(abs(beta.i$lam[beta.i$method == "LASSO"]),2)
    rlasso <- sum(obj$cluster == which.max(obj$centers))/ntotal

    rbsvr <- sum(beta.i[beta.i$method == "BVSR",]$prob > 0.5)/ntotal
    
    return(data.frame(Recover=c(rtpb,rlasso, rbsvr),method = c("TPBayes","LASSO","BVSR"), n.effects=n.effects, std.err=std.err, seed=i.seed))

}))

recover$n.effects <- paste0("n.effects = ",recover$n.effects)
recover$n.effects <- factor(recover$n.effects,levels=unique(recover$n.effects))
recover$std.err <- paste0("sd.err = ",recover$std.err)
recover$std.err <- factor(recover$std.err, levels=unique(recover$std.err ))

recover$method <- factor(recover$method, levels=unique(recover$method ))

p <- ggplot(recover, aes(y=Recover,x=method)) + geom_boxplot(width=0.5, outlier.shape = 20, outlier.size=1) + ylim(0,1) +
     facet_grid(std.err ~ n.effects) + theme(axis.text.x = element_text(angle = 45, hjust = 1))






















file.csv <- paste0("tables/betas_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_to_others.csv")
#write.csv(as.data.frame(beta),file.csv,row.names=F)

p <- ggplot(beta.all[beta.all$std.err==1, ], aes(lam,beta)) + geom_point(shape=1) + facet_grid(method ~ n.effects) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p <- ggplot(beta.all[beta.all$std.err==1, ], aes(lam,beta)) + geom_point(shape=1) + facet_grid(method ~ n.effects) + theme(axis.text.x = element_text(angle = 45, hjust = 1))


################ only tpb test a and b
a=0.1; b=10
beta.all <- c()
pred.all <- c()
  res <- run_sim(param.config[param.config$method=="TPB" & param.config$n.effects == 100,])
    res <- extract_res(res)
    beta.all <- rbind(beta.all, res$betas)
    pred.all <- rbind(pred.all, res$pred)

p <- ggplot(beta.all, aes(lam,beta)) + geom_point(shape=1) + facet_grid(std.err ~ n.effects) + theme(axis.text.x = element_text(angle = 45, hjust = 1))


i=3
        n.effects <- param.config[i,"n.effects"]
        std.err <- param.config[i,"std.err"]
        i.seed <- param.config[i,"seed"]

        set.seed(123 * i.seed)
        
        param <- data.frame(n.effects=n.effects, std.err=std.err,seed=i.seed)

        data=gen_TPBayes_data(std.err = std.err, p=ng, nf=nf, n.effects=n.effects)

results <- TPBayesR(y=data$y[,1:400],x=data$x[,1:400], a=a,b=b,c=c,d=d,g=g,h=h,rsd=123,itr=101, tol = 1e-3, out_dir="results", out_itr = 20)

plot(results$lam,data$beta)
plot(results$z[1,])
sum(which(results$z[1,] > mean(c(min(results$z[1,]),max(results$z[1,])))) %in% which(data$beta != 0))
sum(which(results$z[1,] > 0.2)%in% which(data$beta != 0))

sum((results$z[1,] > mean(c(min(results$z[1,]),max(results$z[1,])))))
ntotal <- sum(data$beta !=0) 
obj <- kmeans(abs(results$lam[1,]),2)
rtpb <- sum(which(obj$cluster == which.max(obj$centers))%in%which(data$beta !=0))/ntotal

which(data$beta !=0)
data$beta[which(data$beta !=0)]

#j = 41
j = 9982
j = 2463
n = 500
x <- cbind(1,data$x[j,])
omega <- solve(t(x) %*% x + diag(1/c(results$sigma_mu, results$phi[j])))
omega0 <- 1/(500+1/results$sigma_mu)
R <- as.numeric(data$y - results$lam %*% data$x + results$lam[,j] * data$x[j,drop=F])
R.bar2 <- mean(R) * mean(R)

ratio = sqrt(det(omega)/omega0) /results$phi[j] * exp(results$psi / 2 *(t(R) %*% x %*% omega %*% t(x) %*% R - omega0*n^2*R.bar2^2))
ratio/(1+ratio)


plot(log(results$phi))

plot(results$z[1,])


results <- run_sim(param.config[1,], method="tpb")
tpb <- extract_res(results)
plot(tpb$beta[,"lam"],tpb$beta[,"beta"])

ntotal <- sum(tpb$beta[,"beta"] !=0 )
obj <- kmeans(abs(tpb$beta[,"lam"]),2)
rtpb <- sum(obj$cluster == which.max(obj$centers))/ntotal

results <- run_sim(param.config[1,], method="SSLASSO")
sslasso <- extract_res(results)
plot(sslasso$beta[,1],sslasso$beta[,3])

file.csv <- paste0("tables/betas_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_to_others.csv")
#write.csv(as.data.frame(beta),file.csv,row.names=F)

p <- ggplot(results.tpb$beta, aes(lam,beta)) + geom_point(shape=1) + facet_grid(std.err ~ n.effects) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot(results.tpb$betas$lam,results.tpb$betas$beta)

results <- run_sim(param.config[1:2,], method="LASSO")
results.tpb <- extract_res(results)

################ comparing tpbayes at a=0.5, b=100 to lasso and bvsr 

##### save betas
#file.csv <- paste0("tables/betas_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_to_others.csv")
#write.csv(as.data.frame(beta),file.csv,row.names=F)

##### save precitions
#file.csv <- paste0("tables/pred_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_to_others.csv")
#write.csv(as.data.frame(pred),file.csv,row.names=F)

##### plot betas
#beta$n.effects <- paste0("n.effects = ",beta$n.effects)
#beta$n.effects <- factor(beta$n.effects,levels=unique(beta$n.effects))

for(i in 1:length(std.err.list)){
    std.err <- std.err.list[i]
    p <- ggplot(beta[beta$std.err ==std.err,], aes(lam,beta)) + geom_point(shape=1) + facet_grid(method ~ n.effects) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

    file.png <- paste0("plots/betas_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_sd_err",std.err,"_to_others.png")

    #pdf(paste("plots/tpbayes_compare_nf100000_a0.5_d100_sd_err_",std.err,".pdf"),width=8,height=6)
    png(file.png,width=1600,height=1200, res = 300)
    print(p)
    dev.off()
}

########## plot predictions
pred$n.effects <- paste0("n.effects = ",pred$n.effects)
pred$n.effects <- factor(pred$n.effects,levels=unique(pred$n.effects))
for(i in 1:length(std.err.list)){
    std.err <- std.err.list[i]
    p <- ggplot(pred[pred$std.err ==std.err,], aes(y.true,y.pred)) + geom_point(shape=1) + facet_grid( method ~ n.effects) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

    file.png <- paste0("plots/predictions_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_sd_err",std.err,"_to_others.png")

    #pdf(paste("plots/tpbayes_compare_nf10000_a0.5_d100_sd_err_",std.err,".pdf"),width=8,height=6)
    png(file.png,width=1600,height=1200, res = 300)
    print(p)
    dev.off()
}

########## plot residuals
pred$residuals <- pred$y.true - pred$y.pred
pred$n.effects <- paste0("n.effects = ",pred$n.effects)
pred$n.effects <- factor(pred$n.effects,levels=unique(pred$n.effects))
pred$std.err <- paste0("sd.err = ",pred$std.err)
pred$std.err <- factor(pred$std.err, levels=unique(pred$std.err ))
p <- ggplot(pred, aes(y=residuals,x=method)) + geom_boxplot(width=0.5, outlier.shape = 20, outlier.size=1) + 
     facet_grid(std.err ~ n.effects) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

file.res <- paste0("plots/residuals_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_to_others.png")

png(file.res,width=1200,height=1200, res = 300)
print(p)
dev.off()

######## count the number of beta that recovered correctly.
file.csv <- paste0("tables/betas_nf",nf,"_comparing_method_",method,"_a",a,"_d",d,"_to_others.csv")
beta <- read.csv(file.csv)

recover <- do.call(rbind,lapply(((1:nrow(param.config))),function(i){
               cat(i,"\n")
  
                n.effects <- param.config[i,"n.effects"]
                std.err <- param.config[i,"std.err"]
                i.seed <- param.config[i,"seed"]
                
                beta.i <- beta[beta$n.effects==n.effects & beta$std.err==std.err & beta$seed==i.seed,]
                
                ntotal <- sum(beta.i[beta.i$method == "TPBayes",]$beta !=0 )

                #rtpb <- sum(beta.i[beta.i$method == "TPBayes",]$prob > 0.99)/ntotal
                obj <- kmeans(abs(beta.i$lam[beta.i$method == "TPBayes"]),2)
                rtpb <- sum(obj$cluster == which.max(obj$centers))/ntotal
                
                obj <- kmeans(abs(beta.i$lam[beta.i$method == "LASSO"]),2)
                rlasso <- sum(obj$cluster == which.max(obj$centers))/ntotal

                rbsvr <- sum(beta.i[beta.i$method == "BVSR",]$prob > 0.5)/ntotal
                
                return(data.frame(Recover=c(rtpb,rlasso, rbsvr),method = c("TPBayes","LASSO","BVSR"), n.effects=n.effects, std.err=std.err, seed=i.seed))

}))

recover$n.effects <- paste0("n.effects = ",recover$n.effects)
recover$n.effects <- factor(recover$n.effects,levels=unique(recover$n.effects))
recover$std.err <- paste0("sd.err = ",recover$std.err)
recover$std.err <- factor(recover$std.err, levels=unique(recover$std.err ))

recover$method <- factor(recover$method, levels=unique(recover$method ))

p <- ggplot(recover, aes(y=Recover,x=method)) + geom_boxplot(width=0.5, outlier.shape = 20, outlier.size=1) + ylim(0,1) +
     facet_grid(std.err ~ n.effects) + theme(axis.text.x = element_text(angle = 45, hjust = 1))




#############################
############################# comparing different parameter settings for tpbayes itself

n.effects.list <- c(20, 50, 100)
std.err.list <- c(1,3,5)

ng <- 1
nf <- 10000

a <- 0.5
b.list <- c(1, 10, 100, 1000, 10000)

results <- do.call(rbind,mclapply((1:length(n.effects.list)),function(i){
        results <- do.call(rbind,mclapply((1:length(std.err.list)),function(k){
            
            results <- do.call(rbind,lapply(1:length(b.list),function(j){
                cat(i,"\n")
                cat(k,"\n")
                set.seed(123)
                n.effects <- n.effects.list[i]
                std.err <- std.err.list[k]
                b <- b.list[j]
                data=gen_TPBayes_data(std.err = std.err, p=ng, nf=nf, n.effects=n.effects)

                results <- TPBayesR(y=data$y[1,1:400],x=data$x[,1:400], a=0.5,b=100,c=0.5,d=b,g=0.5,h=100,rsd=123,itr=101, tol = 1e-3, out_dir="results", out_itr = 20)
                results.tpbayes <- data.frame(lam=as.numeric(results$lam),beta=as.numeric(data$beta[1,]),n.effects=n.effects, std.err=std.err,method="TPBayes", a=a, b=b)
                return(results.tpbayes)
            }))
            return(results)
        },mc.cores = 3))
        return(results)
},mc.cores = 3))

write.csv(as.data.frame(results),"tables/results_nf10000_comparing_tpbayes_tau_to_theta_a0.5_d_vary.csv",row.names=F)

results <- read.csv("tables/results_nf10000_comparing_tpbayes_tau_to_theta_a0.5_d_vary.csv")
#results.bak <- results
#results <- results.bak
results <- as.data.frame(results)
results$n.effects <- paste0("n.effects=",results$n.effects)
results$n.effects <- factor(results$n.effects,levels=unique(results$n.effects))
results$b <- paste0("b=",results$b)

for(i in 1:length(std.err.list)){
    std.err <- std.err.list[i]
    p <- ggplot(results[results$std.err ==std.err,], aes(lam,beta)) + geom_point(shape=1) + facet_grid(b ~ n.effects)

    #pdf(paste("plots/tpbayes_compare_nf10000_a0.5_b100_sd_err_",std.err,".pdf"),width=8,height=6)
    png(paste("plots/tpbayes_compare_nf10000_a0.5_d_vary_sd_err_",std.err,".png"),width=2000,height=1600, res = 300)
    print(p)
    dev.off()
}



}
















