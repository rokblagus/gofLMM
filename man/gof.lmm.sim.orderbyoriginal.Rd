% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_gofLMM.R
\name{gof.lmm.sim.orderbyoriginal}
\alias{gof.lmm.sim.orderbyoriginal}
\title{Goodness-of fit test for LMM, all faster, ordering the residuals by the original fitted values also for the SF and permutation approach}
\usage{
gof.lmm.sim.orderbyoriginal(
  fit,
  std.type = c(1, 2),
  use.correction.for.imbalance = FALSE,
  M = 100,
  verbose = FALSE
)
}
\arguments{
\item{fit}{The result of a call to \code{"nlme"}. The model must be fitted with \code{control=lmeControl( returnObject = TRUE)} and \code{keep.data=TRUE}; ID variable must be numeric and ordered from 1:N !.}

\item{std.type}{Type of standardization to be used for the residuals when constructing the process.
Currently implemeneted options are \code{1} and \code{2} for $S_i=\hat\sigma^{-1/2}I_{n_i}$ and $S_i=\hat{V}_i^{-1/2}$.}

\item{use.correction.for.imbalance}{Logical. use $n_i^{-1/2} S_i$ when standardizing the residuals. Defaults to \code{FALSE}.}

\item{M}{Number of random simulations/sign-flipps/permutations. Defaults to \code{100}.}
}
\value{
KS and CvM pvalues for Pan, Simulation, sign-flip and permutations.
}
\description{
Goodness-of fit test based on cumulative sum stochastic process. Used for simulations. Returns only KS and CvM p-values for all 4 methods and individual as well as cluster specific residuals. Fs not implemented here.
}
\details{
none.
}
\examples{
# simulate some data:
N=50
set.seed(1)
n<-floor(runif(N,min=1,max=15)) #imbalanced
betas<-c(1,1,1,15) #don't change! #the last one is only used whe omit.important.predictor=TRUE
norm.eps<-FALSE
shape=0.5
scale=1
norm.re.intercept<-FALSE
shape.re.intercept=0.5
scale.re.intercept=1
norm.re.slope<-FALSE
shape.re.slope=0.5
scale.re.slope=1
sim.re.slope=FALSE
over.parameterized.model=FALSE #i.e. fit a variable which is not used when generating the data
omit.important.predictor=FALSE
yy<-NA
x22<-NA
id<-NA
x1<-NA
for (gg in 1:N){

  id<-c(id,rep(gg,each=n[gg]))
  x11<-rep(rbinom(1,size=1,prob=0.4),each=n[gg])
  x1<-c(x1,x11)

  if (norm.re.intercept==TRUE) re.int<-rnorm(1,sd=sqrt(2)) else re.int<-rgamma(1,shape=shape.re.intercept,scale=scale.re.intercept)-shape.re.intercept*scale.re.intercept

  b<-rep(re.int,each=n[gg])

  if (norm.re.slope==TRUE) re.slope<-rnorm(1,sd=sqrt(1)) else re.slope<-rgamma(1,shape=shape.re.slope,scale=scale.re.slope)-shape.re.slope*scale.re.slope

  b2<-rep(re.slope,each=n[gg])
  x2<-1:n[gg]
  x4<-runif(n[gg])

  if (norm.eps==TRUE) eps<-rnorm(n[gg]) else eps<-rgamma(n[gg],shape=shape,scale=scale)-shape*scale

  if (sim.re.slope==TRUE) {
    if (omit.important.predictor==FALSE) y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+b2*x2+eps else y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+b2*x2+eps+betas[4]*x4
  } else {
    if (omit.important.predictor==FALSE) y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+eps else y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+eps+betas[4]*x4
  }
  yy<-c(yy,y)
 x22<-c(x22,x2)
}
yy<-yy[-1]
x22<-x22[-1]
x1<-x1[-1]
id<-id[-1]
x4<-runif(sum(n))
aids.art<-data.frame(ptnt=id,outcome=yy,x1=x1,x2=x22,x4=x4)
library(nlme)
fit<-lme(fixed=outcome~ x2+x1:x2, data=aids.art, random=~x2|ptnt,control=lmeControl( returnObject = TRUE),method="REML" )
gof.lmm.sim.orderbyoriginal(fit,std.type=2,use.correction.for.imbalance=FALSE,M=25,verbose=TRUE)
}
\seealso{
\code{\link{gof.lmm.pan}} and \code{\link{gof.lmm}}
}
\author{
Rok Blagus, \email{rok.blagus@mf.uni-lj.si}
}
