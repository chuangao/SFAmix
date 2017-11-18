
SFAmix <- function(Y_TMP_param,nrow_param, ncol_param, a_param,b_param, nf_param, itr_param, LAM_out, EX_out, Z_out, EXX_out, nf_out, out_itr, out_dir,itr_final){
    Y_TMP_param <- as.numeric(as.character(Y_TMP_param))   
    LAM_out <- rep(0,nrow_param*nf_param)
    EX_out <- rep(0,nf_param*ncol_param)
    EXX_out <- rep(0,nf_param*nf_param)
    Z_out <- rep(0,nf_param)
    nf_out <- rep(0,1)
    itr_final <- rep(0,1)
    
    result <- .C ("SFAmix",
                  as.double(Y_TMP_param),as.integer(nrow_param), as.integer(ncol_param), as.double(a_param),as.double(b_param), as.integer(nf_param), as.integer(itr_param), LAM=as.double(LAM_out), EX=as.double(EX_out), Z=as.double(Z_out), EXX=as.double(EXX_out), nf=as.integer(nf_out), as.integer(out_itr), as.character(out_dir),itr_final=as.integer(itr_final))

    nf <- result[['nf']][1]
    itr <- result[['itr_final']][1]
    
    LAM=result[['LAM']]
    LAM <- LAM[1:(nrow_param*nf)]
    LAM <- matrix(LAM,nrow=nrow_param,ncol=nf)

    EX=result[['EX']]
    EX <- EX[1:(ncol_param*nf)]
    EX <- matrix(EX,nrow=nf,ncol=ncol_param)
    
    EXX=result[['EXX']]
    EXX <- EXX[1:(nf*nf)]
    EXX <- matrix(EXX,nrow=nf,ncol=nf)
    
    Z <- 1- result[["Z"]][1:nf]
    
    #return(result)
    return(list(lam=LAM,ex=EX,z=Z,exx=EXX,nf=nf,itr=itr))
    ##return(list(LAM=result$LAM_out,EX=result$EX_out))
}

#' An algorithm for decomposing a high dimensional matrix into the product of a sparse loading matrix, and a dense factor matrix.

#' @author Chuan Gao <chuan.gao.cornell@@gmail.com>

#' @param y matrix to be decmoposed, no missing values are allowed
#' @param nf the number of factors for the algorithm to start with, will be shrank to a smaller number reflecting the number of factors needed to explain the variance, default to 50
#' @param a paramater one for the three parameter beta distribution, default to 0.5 to recapitulate horseshoe
#' @param b paramater two for the three parameter beta distribution, default to 0.5 to recapitulate horseshoe
#' @param itr The maximum number of iterations the algorithm is allowed to run, default to 500
#' @param out_itr (Optional) Iteration number out_itr, the algorithm will write temporary results into the specified directory (see below) every out_itr number of iterations.
#' @param out_dir (Optional) Directory where the algorithm will write temporary results into at the specified iteration number(see above)


#' @return lam: the sparse loading matrix
#' @return ex: the factor matrix
#' @return z: a vector indicating whether the corresponding loading is sparse (value of 1)
#' @return nf: the number of factors learned by the model
#' @return exx: the expected value of the covarance matrix, E(XX^T)
#' @return itr: the number of iterations for the algorithm to converge

#' @examples
#' library(SFAmix)
#' ## simulate data
#' data = gen_SFAmix_data(std=2)

#' ## run algorithm on the simulated data
#' result = SFAmixR(data$y,nf=50,a=0.5,b=0.5,itr=1000)

#' ## calculate a correlation matrix of the estimated loading matrix 
#' ## and the true loading matrix. Ideally, there should be one and 
#' ## only one big correlation value for a given row and column of the 
#' ## correlation matrix
#' cor.est.real = cor(result$lam[,result$z==1],data$lams)

#' ## visulize the correlation matrix
#' image(cor.est.real)

#' @references \url{https://arxiv.org/abs/1310.4792}

SFAmixR <- function(y=y,nf=50,a=0.5,b=0.5,itr=500,out_itr=20,out_dir=NULL){
    out_dir2 = out_dir
    out_dir2 = gsub("/","%",out_dir2)
    if(missing(y)){
        stop("Please check our documentation for the correct usage of this methods, you are missing an input matrix!")
    }
    if(is.null(out_dir)){
       out_dir2 = "NULL"
    }else{
        if(!file.exists(out_dir)){
            stop("Your specified output directory is not found!")
        }
    }
    out_dir2 <- gsub("/","",out_dir2)
    sn = nrow(y)
    dy = ncol(y)
    
    LAM_out <- c()
    EX_out <- c()
    EXX_out <- c()
    Z_out <- c()
    nf_out <- c()
    itr_final <- c()

    result <- SFAmix(y,sn,dy,a,b,nf,itr,LAM_out,EX_out,Z_out, EXX_out, nf_out, out_itr, out_dir2,itr_final)
    
    
    return(result)
}



#' Simulate matrix with dimension of 500 x 200, number of factors is set to 15, where 10 of them being sparse. The sparse loading matrix cotains mostly zeros, and random blocks of nonzero values generated from N(0,std). The dense loading matrix is generated from N(0,std), the factor matrix and the error matrix are generated from N(0,1).

#' @param std standard deviation for the normal distribution

#' @return a list containing the following
#' @return lams: the sparse loading matrix
#' @return lamd: the dense loading matrix
#' @return ex: the factors matrix
#' @return y: the y matrix calculated as y = lam * ex + err
 
gen_SFAmix_data <- function(std=1){
    nf.s <- 10
    nf.d <- 5

    ng <- 500
    ns <- 200

    nf <- nf.s+nf.d


    lams <- matrix(0,nrow=ng,ncol=nf.s)

    lamd <- matrix(rnorm(ng*nf.d,0,std),nrow=ng,ncol=nf.d)

    ex <- matrix(rnorm(nf*ns,0,1),nrow=nf)

    block <- ng/nf.s
    for(i in 1:nf.s){
        ne <- sample(20:40,1)
        lams[(i-1)*block + sample(1:block,ne,replace=F),i] = rnorm(ne,0,std)
    }

    err <- matrix(rnorm(ng*ns),nrow=ng,ncol=ns)

    lam <- as.matrix(cbind(lams,lamd))

    lam <- lam[,sample(1:nf,nf,replace=F)]

    y <- lam %*% ex + err
    
    return(list(y=y,lams=lams,lamd=lamd,ex=ex))
}
