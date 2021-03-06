% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BicMixWrapper.R
\name{BicMixR}
\alias{BicMixR}
\title{An algorithm for decomposing a high dimensional matrix into the product of a sparse loading matrix, and a sparse factor matrix.}
\usage{
BicMixR(y = y, nf = 100, a = 0.5, b = 0.5, itr = 5001,
  rsd = NULL, out_itr = 500, out_dir = NULL, x_method = NULL,
  tol = NULL, qnorm = TRUE)
}
\arguments{
\item{y}{matrix to be decomposed, no missing values are allowed}

\item{nf}{the number of factors for the algorithm to start with, will be shrank to a smaller number reflecting the number of factors needed to explain the variance, default to 100}

\item{a}{paramater one for the three parameter beta distribution, default to 0.5 to recapitulate horseshoe}

\item{b}{paramater two for the three parameter beta distribution, default to 0.5 to recapitulate horseshoe}

\item{itr}{The maximum number of iterations the algorithm is allowed to run, default to 5000}

\item{rsd}{random seed for initializing the parameter values, default to be randomly drawn}

\item{out_itr}{Iteration number out_itr, the algorithm will write temporary results into the specified directory (see below) every out_itr number of iterations. default to 500}

\item{out_dir}{Directory where the algorithm will write temporary results into at the specified iteration number(see above)}

\item{x_method}{whether induce sparsity on the X matrix, take values either "sparse" or "dense". default to "sparse"}

\item{tol}{tolerance threshold for convergence, default to 1e-5}

\item{qnorm}{whether to qq-normalize the gene expression matrix, default to TRUE}
}
\value{
lam: the sparse loading matrix

ex: the factor matrix

z: a vector indicating whether the corresponding loading is sparse (value of 1)

o: a vector indicating whether the corresponding factor is sparse (value of 1)

nf: the number of factors learned by the model

exx: the expected value of the covariance matrix, E(XX^T)
}
\description{
An algorithm for decomposing a high dimensional matrix into the product of a sparse loading matrix, and a sparse factor matrix.
}
\examples{
library(BicMix)
## simulate data, the parameter std specifies the standard error of non-zero entries in the
## loading and factor matrices, where a normal distribution of mean zero
## is assumed for these values.
data = gen_BicMix_data(std=2)
## Visualize the loading matrix
image(t(data$lam),x=1:ncol(data$lam),y=1:nrow(data$lam),xlab="Loadings",ylab="Samples")
## Visualize the factor matrix
image(t(data$ex),x=1:ncol(data$ex),y=1:nrow(data$ex),xlab="Samples",ylab="Factors")
## run algorithm on the simulated data
system("mkdir results")
result = BicMixR(data$y,nf=100,a=0.5,b=0.5,itr=5000,out_dir="results",tol=1e-5,x_method="sparse",rsd=123)
## calculate a correlation matrix of the estimated loading matrix 
## and the true loading matrix. Ideally, there should be one and 
## only one big correlation value for a given row and column of the 
## correlation matrix if the recovered sparse loadings and the true sparse loadings
cor.est.real = cor(result$lam[,result$z==1],data$lams)
## visualize the correlation matrix
image(cor.est.real,x=1:nrow(cor.est.real),y=1:ncol(cor.est.real),
xlab="Recovered loadings",ylab="True loadings")
}
\references{
\url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004791}
}
\author{
Chuan Gao <chuan.gao.cornell@gmail.com>
}
