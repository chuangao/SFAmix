# SFAmix

\name{gen_SFAmix_data}
\alias{gen_SFAmix_data}
\title{Simulate matrix with dimension of 500 x 200, number of factors is set to 15, where 10 of them being sparse. The sparse loading matrix cotains mostly zeros, and random blocks of nonzero values generated from N(0,std). The dense loading matrix is generated from N(0,std), the factor matrix and the error matrix are generated from N(0,1).}
\usage{
gen_SFAmix_data(std = 1)
}
\arguments{
\item{std}{standard deviation for the normal distribution}
}
\value{
a list containing the following

lams: the sparse loading matrix

lamd: the dense loading matrix

ex: the factors matrix

y: the y matrix calculated as y = lam * ex + err
}
\description{
Simulate matrix with dimension of 500 x 200, number of factors is set to 15, where 10 of them being sparse. The sparse loading matrix cotains mostly zeros, and random blocks of nonzero values generated from N(0,std). The dense loading matrix is generated from N(0,std), the factor matrix and the error matrix are generated from N(0,1).
}


\name{SFAmixR}
\alias{SFAmixR}
\title{An algorithm for decomposing a high dimensional matrix into the product of a sparse loading matrix, and a dense factor matrix.}
\usage{
SFAmixR(y = y, nf = 50, a = 0.5, b = 0.5, itr = 500)
}
\arguments{
\item{y}{matrix to be decmoposed, no missing values are allowed}

\item{nf}{the number of factors for the algorithm to start with, will be shrank to a smaller number reflecting the number of factors needed to explain the variance, default to 50}

\item{a}{paramater one for the three parameter beta distribution, default to 0.5 to recapitulate horseshoe}

\item{b}{paramater two for the three parameter beta distribution, default to 0.5 to recapitulate horseshoe}

\item{itr}{The maximum number of iterations the algorithm is allowed to run, default to 500}
}
\value{
lam: the sparse loading matrix

ex: the factor matrix

z: a vector indicating whether the corresponding loading is sparse (value of 1)

nf: the number of factors learned by the model

exx: the expected value of the covarance matrix, E(XX^T)
}
\description{
An algorithm for decomposing a high dimensional matrix into the product of a sparse loading matrix, and a dense factor matrix.
}
\examples{
library(SFAmix)
## simulate data
data = gen_SFAmix_data(std=2)
## run algorithm on the simulated data
result = SFAmixR(data$y,nf=50,a=0.5,b=0.5,itr=1000)
## calculate a correlation matrix of the estimated loading matrix 
## and the true loading matrix. Ideally, there should be one and 
## only one big correlation value for a given row and column of the 
## correlation matrix
cor.est.real = cor(result$lam[,result$z==1],data$lams)
## visulize the correlation matrix
image(cor.est.real)
}
\references{
\url{https://arxiv.org/abs/1310.4792}
}
\author{
Chuan Gao <chuan.gao.cornell@gmail.com>
}
