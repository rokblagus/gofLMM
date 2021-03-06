% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_gofLMM.R
\name{gof.lmm}
\alias{gof.lmm}
\title{Goodness-of fit test for LMM}
\usage{
gof.lmm(
  fit,
  residuals = c("individual", "cluster"),
  std.type = c(1, 2),
  use.correction.for.imbalance = FALSE,
  subset.fix = NULL,
  type = c("simulation", "sign.flip", "permutation"),
  M = 100,
  order.by.original = TRUE,
  force.permutation.with.O = FALSE,
  verbose = FALSE,
  flip.cluster = TRUE,
  use.normal = FALSE,
  use.mammen = FALSE,
  use.sigmoid = FALSE,
  lambda = 0.5,
  transform = TRUE
)
}
\arguments{
\item{fit}{The result of a call to \code{"nlme"}. The model must be fitted with \code{control=lmeControl( returnObject = TRUE)} and \code{keep.data=TRUE}. An error message is returned otherwise. ID variable must be numeric and ordered from 1:N ! Canno't use transofrmations of the outcome variable directly in the formula i.e. lme(sqrt(y)~x) will return p=1!}

\item{residuals}{Residuals to be used when constructing the process. Possible values are \code{"individual"} and \code{"cluster"} for \textit{individual} and \textit{cluster-speciffic} residuals, respectively.}

\item{std.type}{Type of standardization to be used for the residuals when constructing the process.
Currently implemeneted options are \code{1} and \code{2} for S_i=\hat\sigma^{-1/2}I_{n_i} and $S_i=\hat{V}_i^{-1/2}$.}

\item{use.correction.for.imbalance}{Logical. use $n_i^{-1/2} S_i$ when standardizing the residuals. Defaults to \code{FALSE}.}

\item{subset.fix}{Two-sided formula. If nonnull, the process $W^{F^s}$ will be constructed using the variables defined on the RHS of the formula. Deafults to \code{NULL} and the process $W^{F^s}$ is not constructed.}

\item{type}{How to obtain the processes $W^m$. Possible values are \code{"simulation"} for the simulation approach, \code{"sign.flip"} for the sign-flipping approach and \code{"permutation"} for the permutation approach. When using \code{type="permutation"}, sign-flipping will be used by default if not specified otherwise by the argument \code{force.permutation.with.O}.}

\item{M}{Number of random simulations/sign-flipps/permutations. Defaults to \code{100}.}

\item{order.by.original}{Logical. Should the residuals in the the processes $W^m$ be ordered by the original fitted values? Defaults to \code{TRUE}.
Makes sense only for \code{type="sign.flip"} and \code{type="permutation"} since when \code{type="simulation"} the ordering is always based on the original predictions.
It is programmed such that J_i is reestimated at each iteration $m$.}

\item{force.permutation.with.O}{Logical. Should the permutations be used also for the O process? Defaults to \code{FALSE}.}

\item{verbose}{Logical. Print the current status of the test. Can slow down the algorithm, but it can make it feel faster. Defaults to \code{FALSE}.}

\item{flip.cluster}{Logical. Should entire cluster be flipped (i.e. should all subjects within the cluster be multiplied with the same random number). Defaults to \code{TRUE}.}

\item{use.normal}{Lolgical. Use normal random variables instead of sign-flip. Defaultes to \code{FALSE}.}

\item{use.mammen}{Logical. Use Mammen's 2 point dostribution instead of sign-flip. Defaultes to \code{FALSE}. Not in use when \code{use.normal=TRUE}}

\item{use.sigmoid}{Logical. Use sigmoid function instead of the indicator, i.e. smooth the process. Defaults to \code{FALSE}.}

\item{lambda}{Smoothing parameter. Not used when \code{use.sigmoid=FALSE}. Defaults to 0.5.}

\item{transform}{Logical. Should the predictions be transformed to 0,1? Defaults to TRUE.}
}
\value{
An object of class \code{"gofLMM"} for which \code{plot} and \code{summary} functions are available.
}
\description{
Goodness-of fit test based on cumulative sum stochastic process
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
fit.gof<-gof.lmm(fit,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,subset.fix=outcome~x2,type= "simulation" ,M=25,order.by.original=FALSE,force.permutation.with=FALSE,verbose=TRUE)
plot.gofLMM(fit.gof,type=2,subset.M=NULL,xlab="",main="Example")
summary.gofLMM(fit.gof)

library(nlme)
data(Orthodont)
Orthodont$Subject<- rep(1:27,each=4)
fm1<-lme(distance~age,random=~1|Subject,data=Orthodont,control=lmeControl( returnObject = TRUE),method="REML")
gof.fm1<-gof.lmm(fm1,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,subset.fix=NULL,type= "sign.flip" ,M=500,order.by.original=FALSE,force.permutation.with.O=FALSE,verbose=TRUE)
plot.gofLMM(gof.fm1,type=2,subset.M=NULL,xlab="",main="Orthodont, model 1")
summary.gofLMM(gof.fm1)

fm1.1<-lme(distance~age,random=~age|Subject,data=Orthodont,control=lmeControl( returnObject = TRUE),method="REML")
gof.fm1.1<-gof.lmm(fm1.1,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,subset.fix=NULL,type= "sign.flip" ,M=500,order.by.original=FALSE,force.permutation.with.O=FALSE,verbose=TRUE)
plot.gofLMM(gof.fm1.1 ,type=2,subset.M=NULL,xlab="",main="Orthodont, model 1.1")
summary.gofLMM(gof.fm1.1)

fm2<-lme(distance~age+Sex,random=~1|Subject,data=Orthodont,control=lmeControl( returnObject = TRUE),method="REML")
gof.fm2<-gof.lmm(fm2,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,subset.fix=distance~age,type= "sign.flip" ,M=500,order.by.original=FALSE,force.permutation.with=FALSE,verbose=TRUE)
plot.gofLMM(gof.fm2,type=2,subset.M=NULL,xlab="",main="Orthodont, model 2")
summary.gofLMM(gof.fm2)

fm2.1<-lme(distance~age*Sex,random=~1|Subject,data=Orthodont,control=lmeControl( returnObject = TRUE),method="REML")
gof.fm2.1<-gof.lmm(fm2.1,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,subset.fix=NULL,type= "sign.flip" ,M=500,order.by.original=FALSE,force.permutation.with.O=FALSE,verbose=TRUE)
plot.gofLMM(gof.fm2.1,type=2,subset.M=NULL,xlab="",main="Orthodont, model 2.1")
summary.gofLMM(gof.fm2.1)
}
\seealso{
\code{\link{gof.lmm.pan}}, \code{\link{plot.gofLMM}} and  \code{\link{summary.gofLMM}}
}
\author{
Rok Blagus, \email{rok.blagus@mf.uni-lj.si}
}
