
cal_score_sparse <- function(sigma, precis = FALSE){

    sigma <- abs(sigma)
    nr <- nrow(sigma)
    nc <- ncol(sigma)

    r1 <- apply(sigma,1,function(x){
        xm <- max(x)
        x <- x[x!=max(x)]
        return(xm - mean(x[x>=mean(x)]))
    })
    r1 <- as.numeric(r1)

    r2 <- 0
    if(nrow(sigma) > 1){
        r2 <- apply(sigma,2,function(x){
            xm <- max(x)
            x <- x[x!=max(x)]
            return(xm - mean(x[x>=mean(x)]))  
        })
    }
    r2 <- as.numeric(r2)

    nf <- min(nc,nr)
    o1 <- order(r1,decreasing=T)    
    o2 <- order(r2,decreasing=T)

    r <- (mean(r1)+mean(r2))/2
    if(precis){
        r <- (mean(r1[o1[1:nf]]) + mean(r2[o2[1:nf]]))/2
    }
    return(r)
}

cal_score_dense <- function(m1,m2,precis=FALSE){
    
    m1 <- apply(m1,2,function(x){scale(x)})
    m2 <- apply(m2,2,function(x){scale(x)})
    if(precis){
        cr <- cor(m1,m2)
        cr <- abs(cr)
        nr <- nrow(cr)
        nc <- ncol(cr)
        nf <- min(nr,nc)
        maxr <- apply(cr,1,max)
        maxc <- apply(cr,2,max)
        o.r <- order(maxr,decreasing=T)
        o.c <- order(maxc,decreasing=T)
        if(nr > nc){
            m1 <- m1[,o.r[1:nc]]
        }else{
            m2 <- m2[,o.c[1:nr]]
        }
    }
    sigma <- m1 %*% t(m1) - m2 %*% t(m2)
    r <- sum(diag(sigma * sigma))/nrow(sigma)/nrow(sigma)
    return(r)
}

#m1 = matrix(rnorm(10*3),nrow=10)
#m2 = matrix(rnorm(10*5),nrow=10)


# gen_TPBayes_data <- function(std.err = 1, std=1, p=1, nf=1000, n.effects = 30, ns=500, shape=1){

#     ng <- p
 
#     nf <- nf

#     n.effects <- n.effects

#     lam <- matrix(rep(0,nf * ng),nrow=ng)
#     lam <- apply(lam, 1, function(x){       
#         index.lam <- sample(1:nf,n.effects,replace=F)
#         x[index.lam] <- rgamma(length(index.lam),shape,1)
#         return(x)
#     })    
#     lam <- t(lam)
    
#     ex <- matrix(rnorm(ns*nf,0,std),nrow=nf,ncol=ns)
 
#     err <- matrix(rnorm(ng*ns,0,std.err),nrow=ng,ncol=ns)

#     y <-  lam %*% ex + err
    
#     return(list(y=y,beta=lam,x=ex))
# }


# gen_bulk_data <- function(param.config, output.dir = "data"){
#     mclapply(1:nrow(param.config),function(i){
#         #cat(i,"\n")
#         #write.table(i,file.path(output.dir,"itr.txt"))
        
#         n.effects <- param.config[i,"n.effects"]
#         std.err <- param.config[i,"std.err"]
#         i.seed <- param.config[i,"seed"]
#         method <- param.config[i,"method"]
#         std.effect <- param.config[i,"std.effect"]
#         b <- param.config[i,"b"]
#         dense <- as.logical(as.character(param.config[i,"dense"]))

#         param <- data.frame(n.effects=n.effects, std.err=std.err,seed=i.seed, std.effect = std.effect,dense=dense)
        
#         data=gen_SFA_data(std=2, rsd = i.seed, std.err=std.err, n.effects = 30, nfs = 20, ng=1000,ns=500)
#         file.name <- paste(paste0(names(param),param[1,]),collapse="_")
        
#         write.table(data$y,file.path(output.dir,paste0(file.name,".txt")),col.names=F,row.names=F)
#  },mc.cores = 10)
# }

##### inputDir is for method like KSVD that require reading data file from hard drive
##### output.Dir is for method like KSVD that require writting results to hard drive.

run_sim <- function(param.config, itr=2001, inputDir=NULL, outputDir = NULL, nfs=20, nf=30, ng= 1000, ns=500, min_fac=20,step=2){
    results<- mclapply(((1:nrow(param.config))),function(i){
        #cat(i,"\n")
        #write.table(i,file.path(output.dir,"itr.txt"))

        inputDir='/Users/cg253/data_BicMix2'
        
        
        nfs <- nfs
        nf <- nf

        nfd <- nf - nfs

        ng <- ng
        ns <- ns
        
        n.effects <- param.config[i,"n.effects"]
        std.err <- param.config[i,"std.err"]
        i.seed <- param.config[i,"seed"]
        method <- as.character(param.config[i,"method"])
        std.effect <- param.config[i,"std.effect"]
        b <- param.config[i,"b"]
        dense <- as.logical(as.character(param.config[i,"dense"]))

        param <- data.frame(n.effects=n.effects, std.err=std.err,seed=i.seed, method = method, std.effect = std.effect, b=b, dense=dense,stringsAsFactors =F)
        
        print(i)
        print(param)

        data=gen_SFA_data(std=2, rsd = i.seed, std.err=std.err, n.effects = nf, nfs=nfs, ng=ng,ns=ns,dense=dense)

        file.name <- paste(paste0(names(param),param[1,]),collapse="_")
        inputFile=paste0(file.name,".txt")
        write.table(data$y,file.path(inputDir,inputFile),col.names=F,row.names=F)

        #param2 <- data.frame(n.effects=n.effects, std.err=std.err,seed=i.seed, std.effect = std.effect, dense=dense)
        
        #file.name <-  paste(paste0(names(param2),param2[1,]),collapse="_")    
     
        c=a;g=a;
        d=round(b/2);h=round(b/2);

        #b = 10000000
        lams <- NULL
        lamd <- NULL
        z <- NULL
        lam <- NULL
        if(method == "BicMix2"){
            results <- c()
            #nf.input <- nfs
            #if(dense){
            #    nf.input = nf
            #}
            results <- BicMixR2(y=data$y, nf = 100, a=a,b=b,c=c,d=d,g=g,h=h,rsd=12345,itr=itr, tol = 0.0001, out_dir=results.path, out_itr = 500,min_fac=min_fac,step=step)

            # results <- BicMixR2(y=data$y, nf = 50, a=a,b=b,c=c,d=d,g=g,h=h,rsd=12345,itr=itr, tol = 0.0001, out_dir=results.path, out_itr = 500)


            psi <- results$psi
            bf <- results$bf
            z <- results$z

            #bf.adj <- bf/sqrt(mean(psi))
          
            #count.prob <- apply(bf,2,function(x){return(sum(x > 1000)/length(x))})     
            count.prob <- apply(z,2,function(x){return(sum(x > 0.5)/length(x))})

            index.sparse <- count.prob < 0.5 & count.prob != 0
            lams <- results$lam[,index.sparse,drop=F]
            lamd <- results$lam[,!index.sparse,drop=F]

            # count.sparse <- apply(results$bf,2,function(x){
            #     obj <- kmeans(abs(x),2)
            #     count <- sum(obj$cluster == which.max(obj$centers))/nrow(results$lam)
            #     return(count)
            # })
            # lams <- results$lam[,count.sparse < 0.1,drop=F]
            # lamd <- results$lam[,count.sparse >= 0.1,drop=F]

        }else if(method == "SPCA"){
            if(dense){
                sp<- SPC(t(data$y),sumabsv=4,K=nf,niter=100)
                lamd <- sp$v[,1:nfd,drop=F]
                lams <- sp$v[,(nfd+1):nf,drop=F]
            }else{
                sp <- SPC(t(data$y),sumabsv=4,K=nfs,niter=100)
                lams <- sp$v[,1:(nfs),drop=F]
            }
            
        }else if(method == "KSVD"){
            outputDir='/Users/cg253/results_BicMix2/KSVD'
            res <- c()
            if(dense){
                res <- ksvd(inputDir=inputDir,inputFileName=file.name,n=ns,p=ng,k=nf,outputDir=outputDir,scriptFile='/Users/cg253/BicMix2/KSVD/mine_sim.m',matlabWhere='/Applications/Matlab_R2019b.app/bin/matlab')
                #inputDir=NULL,inputFileName=NULL,n=NA,p=NA,k=NA, outputDir=NULL,scriptFile=NULL,matlabWhere
            }else{
                res <- ksvd(inputDir=inputDir,inputFileName=file.name,n=ns,p=ng,k=nfs,outputDir=outputDir,scriptFile='/Users/cg253/BicMix2/KSVD/mine_sim.m',matlabWhere='/Applications/Matlab_R2019b.app/bin/matlab')
            }
            
            lam <- res$lam
            index.lam <- order(apply(lam,2,var))

            lams <- lam[,index.lam[1:(nfs)],drop=F]
            if(dense){  
                lamd <- lam[,index.lam[(nfs+1):nf],drop=F]
            }
            
        }else if(method == "IFA"){ 
            outputDir='/Users/cg253/results_BicMix2/IFA'
            res <- c()
            if(dense){
                res <- ifa(inputDir=inputDir,inputFileName=file.name,n=ns,p=ng,k=nf,outputDir=outputDir,scriptFile='/Users/cg253/BicMix2/IFA/IFA_chuan.m',matlabWhere='/Applications/Matlab_R2019b.app/bin/matlab')
            }else{
                res <- ifa(inputDir=inputDir,inputFileName=file.name,n=ns,p=ng,k=nfs,outputDir=outputDir,scriptFile='/Users/cg253/BicMix2/IFA/IFA_chuan.m',matlabWhere='/Applications/Matlab_R2019b.app/bin/matlab')
            }
            
            lam <- res$lam
            lam <- lam[,ncol(lam):1,drop=F]
            count <- apply(lam,2,function(x){return(sum(x!=0))})
            lam <- lam[,count>0,drop=F]
            ncol <- ncol(lam)
            if(ncol <= nfs){
                lams <- lam[,1:ncol,drop=F]
            }else{
                lams <- lam[,1:nfs,drop=F]
            }
            if(dense){  
                if(ncol > nfs){
                    lamd <- lam[,(nfs+1):ncol,drop=FALSE]
                }else{
                    lamd <- NULL
                }
            }         
        }else if(method == "BFRM"){  
            outputDir='/Users/cg253/results_BicMix2/BFRM'
            #system(paste0("mkdir -p ",outputDir))
            if(dense){
                res <- bfrm(inputDir=inputDir, inputFileName=file.name,n=ns,p=ng,k=nf, outputDir=outputDir,scriptFile='/Users/cg253/BicMix2/BFRM/parameters.txt',bfrmWhere='/Users/cg253/BicMix2/BFRM/bfrm')
            }else{
                res <- bfrm(inputDir=inputDir, inputFileName=file.name,n=ns,p=ng,k=nfs, outputDir=outputDir,scriptFile='/Users/cg253/BicMix2/BFRM/parameters.txt',bfrmWhere='/Users/cg253/BicMix2/BFRM/bfrm')
            }
            
            lam <- res$lam
            index.lam <- order(apply(lam,2,var))
    
            lams <- lam[,index.lam[1:(nfs)],drop=F]
            if(dense){  
                lamd <- lam[,index.lam[(nfs+1):nf],drop=F]
            }

        }

        else if(method == "BicMix"){
            results <- BicMixR(y=data$y, nf = 50, out_dir=results.path, rsd = 123, tol=1e-10, itr = 1001)
            lam <- results$lam
            z <- results$z
            lams <- lam[,z>0.8,drop=F]
            lamd <- lam[,z<=0.8,drop=F]
        }

        nfs.o <- ifelse(is.null(lams),NA,ncol(lams))
        nfd.o <- ifelse(is.null(lamd),NA,ncol(lamd))

        score.sparse <- NA
        score.sparse.precis <- NA

        if(length(lams) == 0){
            score.sparse <- NA
        }else{
            corr.sparse <- cor(lams,data$lams)
            score.sparse <- cal_score_sparse(corr.sparse)
            score.sparse.precis <- cal_score_sparse(corr.sparse,precis=TRUE)
        }
        score.sparse.all <- data.frame(param,score.sparse = score.sparse,score.sparse.precis = score.sparse.precis, nfs=nfs.o,nfd=nfd.o)

        score.dense <- NA
        score.dense.precis <- NA

        if(dense){
            if(length(lamd) == 0){
                score.dense <- NA
            }else{
                score.dense <- cal_score_dense(lamd,data$lamd)
                score.dense.precis <- cal_score_dense(lamd,data$lamd, precis=TRUE)

            }
        }
        score.dense.all <- data.frame(param,score.dense = score.dense,score.dense.precis = score.dense.precis, nfs=nfs.o,nfd=nfd.o)
       
        if(dense){
            return(list(score.sparse=score.sparse.all, score.dense=score.dense.all))
        }else{
            return(list(score.sparse=score.sparse.all, score.dense=data.frame(param,score.dense = NA,score.dense.precis = NA, nfs=nfs.o,nfd=nfd.o)))
        }
        
        
    },mc.cores = 10)
    return(results)
}

extract_res <- function(results){
    
    ##### extract betas
    score.sparse <- do.call(rbind,lapply(results,function(x){
        return(x$score.sparse)     
    }))

    ##### extract predicted values
    score.dense <- do.call(rbind,lapply(results,function(x){
        return(x$score.dense)     
    }))

    # lams.recover <- do.call(rbind,lapply(results,function(x){
    #     return(x$lams.recover)     
    # }))

    # lams.true <- do.call(rbind,lapply(results,function(x){
    #     return(x$lams.true)     
    # }))

    # return(list(score.sparse = score.sparse, score.dense = score.dense, lams.recover = lams.recover, lams.true = lams.true))
    return(list(score.sparse = score.sparse, score.dense = score.dense))
}

cal_roc <- function(df){
    obj <- roc(df$prob.true,df$prob)
    return(data.frame(fp=1-obj$specificities,tp=obj$sensitivities))
}

cal_quant <- function(df){
    #df <- df[order(df$fp),]
    #print(df[1,])
    prob=seq(0,1,0.01)
    fp.map.to.prob <- cut(df$fp,prob,include.lowest=TRUE)
    index.non.dup <- !duplicated(fp.map.to.prob)
    tp.map.to.prob <- df$tp[index.non.dup]
    fp.tmp = fp.map.to.prob[index.non.dup]

    d = data.frame(
        fp = prob[fp.tmp],
        tp = tp.map.to.prob
    )
    return(d)
}



gen_SFAmix_data <- function(std=2, rsd = 123, std.err=1, n.effects = 30, nfs = 20, ng = 1000, ns=200, dense=T){
    set.seed(rsd)
    
    nf <- n.effects
    nfd <- nf - nfs

    lam <- c()
    lams <- c()
    lamd <- c()

    lams <- matrix(0,nrow=ng,ncol=nfs)

    ########## simulate lam
    block <- 10
    for(i in 1:nfs){
        start <- sample(1:(ng-block),1)
        index <- start + sample(1:block,block,replace=T)
        lams[index,i] = rnorm(block,1,std)
        lams[sample(index,block/2),i] = rnorm(block/2,-1,std)
    }
    
    if(dense){
        l <- ng*nfd
        lamd.tmp <- rnorm(l,1,std)
        lamd.tmp[sample(1:l,l/2)] <- rnorm(l/2,-1,std)
        lamd <- matrix(lamd.tmp,nrow=ng,ncol=nfd)
        lam <- cbind(lams,lamd)
        lam <- lam[,sample(1:nf,nf,replace=F)]
    }else{
        lam <- lams
        lam <- lam[,sample(1:nfs,nfs,replace=F)]
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
    p <- ggplot(m,aes(x=Var1,y=Var2)) + geom_tile(aes(fill=value)) + coord_flip() +theme_minimal() + theme(legend.position="none",axis.title = element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),plot.title = element_text(face = "bold", hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(title) + scale_fill_gradient2(low="blue",mid="white",high="red") 
    #scale_fill_gradient2(low="green",mid="white",high="red") 
    #panel.border = element_rect(colour = "black", fill=NA, size=1)
}

sfa_scheme <- function(){
    #data=gen_SFA_data2(std=2, rsd = 123, std.err=1, n.effects = 5, nfs=3, ng=10,ns=5,dense=TRUE)
    data=gen_SFA_data3(std=0.5, rsd = 123, std.err=1, n.effects = 15, nfs=10, ng=50,ns=20,dense=TRUE)

    plam <- mine_heatmap(data$lam,"Loading")
    pex <- mine_heatmap(data$ex,"Factor")
    py <- mine_heatmap(data$y,"Y")
    perr <- mine_heatmap(data$err,"Error")
    return(list(plam = plam, pex = pex, py = py, perr = perr))
}
