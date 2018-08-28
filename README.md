## SFAmix

###### SFAmix is a sparse matrix decomposition tool. Given a matrix Y with dimension of P by N, SFAmix decompose it into the product of a sparse matrix LAM and a dense matrix X

###### This is the C++ implementation of SFAmixC wrapped in R. It is compiled single threaded. If you want to run multiple threaded, please check the SFAmixC package 

## Use devtools to install in R
`library(devtools)` <br/>
`install_github("chuangao/SFAmix")` <br/>

## Install from source
If install_github command fails (I found that install_github can't resolve the namespace that I specificy the .C all), then clone library into one of your local directory, then install from source <br/>

`git clone https://github.com/chuangao/SFAmix` <br/>
`R CMD INSTALL SFAmix` <br/>

## Usage

SFAmixR(y = y, nf = 100, itr = 5000) <br/>

**Please no headers in the input matrix, no missing values, just pure numbers, ideally quantile normalized** <br/>
**Also no corrections of confounding beforehand, SFAmix will handle that in the dense components** <br/>
**For a gene expression matrix, it is prefered that each gene is a row and each sample is a column** <br/> 

### Arguments
**y** matrix to be decmoposed, no missing values are allowed <br/>
**nf** the number of factors for the algorithm to start with, will be shrank to a smaller number reflecting the number of factors needed to explain the variance, default to 50 <br/>
**itr** The maximum number of iterations the algorithm is allowed to run, default to 5000 <br/>
**out_itr** (Optional) the algorithm will write temporary results into the specified directory (see below) every out_itr number of iterations <br/>
**out_dir** (Optional) Directory where the algorithm will write temporary results into at the specified iteration number(see above) <br/>

### Value
**lam** the sparse loading matrix <br/>
**ex** the factor matrix <br/>
**z** a vector indicating whether the corresponding loading is sparse (value of 1) <br/>
**o** a vector indicating whether the corresponding factor is sparse (value of 1) <br/>
**nf** the number of factors learned by the model <br/>
**exx** the expected value of the covarance matrix, E(XX^T) <br/>
**itr** the number of iterations for the algorithm to converge <br/>

### Examples
library(SFAmix) <br/>
\# simulate data <br/>
data = gen_SFAmix_data(std=2) <br/>
\# run algorithm on the simulated  <br/>
result = SFAmixR(data$y,nf=50,itr=5000) <br/>
\# calculate a correlation matrix of the estimated loading matrix <br/>
\# and the true loading matrix. Ideally, there should be one and only one big correlation value for a given row and column of the correlation matrix <br/>
cor.est.real = cor(result$lam[,result$z==1],data$lams) <br/>
\# visulize the correlation matrix <br/>
image(cor.est.real) <br/>

### Documentation
Please refer to SFAmix.pdf for more usage details <br/>

### References
This work is detailed in <br/>
A latent factor model with a mixture of sparse and dense factors to model gene expression data with confounding effects<br/>
https://arxiv.org/abs/1310.4792
