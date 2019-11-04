


##############################

##aux functions

#' Internal function
#' @keywords internal


"%^%" <- function(x, n)   with(eigen(x), vectors %*% (values^n * t(vectors)))



#' Internal function
#' @keywords internal
#' @param ... additional arguments passed to from or to other methods

my.plot.gofLMM.part<-function(W,Wm,type=c(1,2),y,ym,...){

if (type==1) {x<-1:length(W);xm<-list(); for (ii in 1:length(Wm)) { xm[[ii]]<-1:length(W) }} else {x<-y[order(y)];xm=lapply(ym,function(x) x[order(x)] )}

if (!is.list(ym)) xm<-x
ylim.min<-min( W,min(unlist(Wm)) )
ylim.max<-max( W,max(unlist(Wm)) )

xlim.min<-min( min(x),min(unlist(xm)) )
xlim.max<-max( max(x),max(unlist(xm)) )

plot(x,W,col="white",type="s",ylim=c(ylim.min,ylim.max),xlim=c(xlim.min,xlim.max),...)
for (ii in 1:length(Wm)){
  if (is.list(ym)) lines(xm[[ii]],Wm[[ii]],type="s",col="lightgray",...) else lines(xm,Wm[[ii]],type="s",col="lightgray",...)
}
lines(x,W,type="s",...)

}



#' Internal function
#' @keywords internal

CvM<-function(x) {sum(x**2)}

#' Internal function
#' @keywords internal

KS<-function(x) {max(abs(x))}

#' Internal function
#' @keywords internal

p.val<-function(testStat,testStatm) {(sum( testStatm>= testStat )+1)/(length(testStatm)+1)}

#' Internal function
#' @keywords internal

test.stat.p.val<-function(W,Wm){
  ks<-KS(W)
  cvm<-CvM(W)

  ksm<-unlist(lapply(Wm,KS))
  cvmm<-unlist(lapply(Wm,CvM))

  ks.p<-p.val(ks,ksm)
  cvm.p<-p.val(cvm,cvmm)

  res<-rbind(c(ks,ks.p),c(cvm,cvm.p))
  colnames(res)<-c("TestStat","p.value")
  rownames(res)<-c("KS","CVM")
  res
}



#' Internal function
#' @keywords internal


get.sim.proc<-function(fit, residuals ,std.type ,use.correction.for.imbalance ,subset.fix ,order.by.original ,or.original.fitted.I ,or.original.fitted.P ,or.original.fitted.S,original.fitted.I ,original.fitted.P ,original.fitted.S,n,N,x,ZZ,id ){


resI<-residuals(fit, level = 1  )


resP<-residuals(fit, level = 0  )

if (order.by.original==TRUE) estI<-original.fitted.I  else estI<-fitted(fit,level=1)
if (order.by.original==TRUE) estP<-original.fitted.P  else estP<-fitted(fit,level=0)

if (order.by.original==TRUE)  orI<-or.original.fitted.I  else orI<-order(estI)
if (order.by.original==TRUE)  orP<-or.original.fitted.P  else orP<-order(estP)


vc<-VarCorr(fit)
sigma.est<-as.numeric(vc[nrow(vc),1])

D<-getVarCov(fit)

beta.f<-fixef(fit)

V<-list()
V.i<-list()
Z<-list()

H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
for (gg in 1:N){
if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
I<-diag(rep(1),n[[gg]])
V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
V.i[[gg]]<-V[[gg]]%^%(-1)
if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
}

H.i<-solve(H)


J<-list()
A<-list()
B<-list()

res.i.c<-resI


for (gg in 1:N){


if (n[gg]!=1) A[[gg]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else A[[gg]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

if (n[gg]!=1) B[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else B[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])



I<-diag(rep(1,n[gg]))

if (residuals=="individual") J[[gg]]<-sigma.est*V.i[[gg]]-(A[[gg]])%*%ginv(B[[gg]])%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]] else J[[gg]]<-I-(A[[gg]]+B[[gg]])%*%ginv(B[[gg]])%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]


if (residuals=="individual") res.i.c[id==gg]<- J[[gg]]%*% resI[id==gg] else  res.i.c[id==gg]<- J[[gg]]%*% resP[id==gg]



}



V.ii.inv<-list()


if (residuals=="individual") res.i.c2<-resI else res.i.c2<-resP


resIst<-NA
resPst<-NA
for (gg in 1:N){
I<-diag(rep(1,n[gg]))

V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)




if (std.type==2) Si<-V.ii.inv[[gg]] else Si<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
if (use.correction.for.imbalance==TRUE) Si<-Si/sqrt(n[gg])

resPMpC<-matrix(res.i.c[id==gg],ncol=1,nrow=n[gg],byrow=F)
resPMpC2<-Si%*%resPMpC
resPMpC2<-resPMpC2

resIst<-c(resIst,resPMpC2)


resPMpCP<-matrix(res.i.c2[id==gg],ncol=1,nrow=n[gg],byrow=F)
resPMpC2P<-Si%*%resPMpCP
resPMpC2P<-resPMpC2P

resPst<-c(resPst,resPMpC2P)

}


resIst<-resIst[-1]
 resPst<-resPst[-1]


resoI2<-resIst[orI]
 t01<- estI

for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
ig<-which(round(t01[orI],10)==round(ii,10))
resoI2[ig]<-sum(resoI2[ig])/length(ig)
}

WI2<-1/sqrt(N )*cumsum(resoI2)

resoP2<-resPst[orP]
 t01P<- estP
for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
ig<-which(round(t01P[orP],10)==round(ii,10))
resoP2[ig]<-sum(resoP2[ig])/length(ig)
}

WP2<-1/sqrt(N )*cumsum(resoP2)


##for Fs:
if (!is.null(original.fitted.S)){

if (order.by.original==FALSE){
x.subset<-model.matrix(subset.fix, data=fit$data   )
cfs.fix.sub<-fixef(fit)[colnames(x.subset)]

estS<-x.subset%*%cfs.fix.sub
orS<-order(estS)
} else {
estS<-original.fitted.S
orS<-or.original.fitted.S
}

resoP22<-resPst[orS]
 t01P<- estS
for (ii in as.numeric(names(table(t01P[orS]))[which(table(t01P[orS])>1)])){
ig<-which(round(t01P[orS],10)==round(ii,10))
resoP22[ig]<-sum(resoP22[ig])/length(ig)
}

WP2s<-1/sqrt(N )*cumsum(resoP22)

list(WI2,WP2,WP2s,estI,estP,estS)
} else list(WI2,WP2,estI,estP)

}




#' Internal function
#' @keywords internal

get.sim.proc.O<-function(fit, residuals ,std.type ,use.correction.for.imbalance ,  order.by.original ,or.original.fitted.I , original.fitted.I  ,n,N,x,ZZ,id ){



resI<-residuals(fit, level = 1  )


resP<-residuals(fit, level = 0  )

if (order.by.original==TRUE)  estI<-original.fitted.I  else estI<-fitted(fit,level=1)


if (order.by.original==TRUE)  orI<-or.original.fitted.I  else orI<-order(estI)



vc<-VarCorr(fit)
sigma.est<-as.numeric(vc[nrow(vc),1])

D<-getVarCov(fit)

beta.f<-fixef(fit)

V<-list()
V.i<-list()
Z<-list()

H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
for (gg in 1:N){
if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
I<-diag(rep(1),n[[gg]])
V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
V.i[[gg]]<-V[[gg]]%^%(-1)
if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
}

H.i<-solve(H)


J<-list()
A<-list()
B<-list()

res.i.c<-resI


for (gg in 1:N){


if (n[gg]!=1) A[[gg]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else A[[gg]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

if (n[gg]!=1) B[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else B[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])



I<-diag(rep(1,n[gg]))

if (residuals=="individual") J[[gg]]<-sigma.est*V.i[[gg]]-(A[[gg]])%*%ginv(B[[gg]])%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]] else J[[gg]]<-I-(A[[gg]]+B[[gg]])%*%ginv(B[[gg]])%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]


if (residuals=="individual") res.i.c[ id==gg]<- J[[gg]]%*% resI[ id==gg] else  res.i.c[ id==gg]<- J[[gg]]%*% resP[ id==gg]



}



V.ii.inv<-list()
V.ii<-list()



resIst<-NA

for (gg in 1:N){
I<-diag(rep(1,n[gg]))

V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
V.ii[[gg]]<-V[[gg]]%^%(0.5)



if (std.type==2) Si<-V.ii.inv[[gg]] else Si<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
if (use.correction.for.imbalance==TRUE) Si<-Si/sqrt(n[gg])

resPMpC<-matrix(res.i.c[id==gg],ncol=1,nrow=n[gg],byrow=F)
resPMpC2<-Si%*%resPMpC
resPMpC2<-resPMpC2

resIst<-c(resIst,resPMpC2)


}


resIst<-resIst[-1]



resoI2<-resIst[orI]
 t01<- estI

for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
ig<-which(round(t01[orI],10)==round(ii,10))
resoI2[ig]<-sum(resoI2[ig])/length(ig)
}

WI2<-1/sqrt(N )*cumsum(resoI2)

list(WI2,estI)

}


#' Internal function
#' @keywords internal


get.sim.proc.F<-function(fit, residuals ,std.type ,use.correction.for.imbalance ,subset.fix ,order.by.original ,  or.original.fitted.P ,or.original.fitted.S, original.fitted.P ,original.fitted.S,n,N,x,ZZ,id ){




resI<-residuals(fit, level = 1  )


resP<-residuals(fit, level = 0  )


if (order.by.original==TRUE)  estP<-original.fitted.P else  estP<-fitted(fit,level=0)


if (order.by.original==TRUE)  orP<-or.original.fitted.P  else orP<-order(estP)


vc<-VarCorr(fit)
sigma.est<-as.numeric(vc[nrow(vc),1])

D<-getVarCov(fit)

beta.f<-fixef(fit)

V<-list()
V.i<-list()
Z<-list()

H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
for (gg in 1:N){
if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
I<-diag(rep(1),n[[gg]])
V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
V.i[[gg]]<-V[[gg]]%^%(-1)
if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
}

H.i<-solve(H)




V.ii.inv<-list()


if (residuals=="individual") res.i.c2<-resI else res.i.c2<-resP



resPst<-NA
for (gg in 1:N){
I<-diag(rep(1,n[gg]))

V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)




if (std.type==2) Si<-V.ii.inv[[gg]] else Si<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
if (use.correction.for.imbalance==TRUE) Si<-Si/sqrt(n[gg])



resPMpCP<-matrix(res.i.c2[id==gg],ncol=1,nrow=n[gg],byrow=F)
resPMpC2P<-Si%*%resPMpCP
resPMpC2P<-resPMpC2P

resPst<-c(resPst,resPMpC2P)

}



 resPst<-resPst[-1]


resoP2<-resPst[orP]
 t01P<- estP
for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
ig<-which(round(t01P[orP],10)==round(ii,10))
resoP2[ig]<-sum(resoP2[ig])/length(ig)
}

WP2<-1/sqrt(N )*cumsum(resoP2)


##for Fs:
if (!is.null(original.fitted.S)){

if (order.by.original==FALSE){
x.subset<-model.matrix(subset.fix, data=fit$data   )
cfs.fix.sub<-fixef(fit)[colnames(x.subset)]

estS<-x.subset%*%cfs.fix.sub
orS<-order(estS)
} else {
estS<-original.fitted.S
orS<-or.original.fitted.S
}

resoP22<-resPst[orS]
 t01P<- estS
for (ii in as.numeric(names(table(t01P[orS]))[which(table(t01P[orS])>1)])){
ig<-which(round(t01P[orS],10)==round(ii,10))
resoP22[ig]<-sum(resoP22[ig])/length(ig)
}

WP2s<-1/sqrt(N )*cumsum(resoP22)

list( WP2,WP2s, estP,estS)
} else list( WP2,estP)

}












#######main function

#' Goodness-of fit test for LMM
#'
#' Goodness-of fit test based on cumulative sum stochastic process

#'
#' @param fit The result of a call to \code{"nlme"}. The model must be fitted with \code{control=lmeControl( returnObject = TRUE)} and \code{keep.data=TRUE}. An error message is returned otherwise. ID variable must be numeric and ordered from 1:N ! Cannot use transofrmations of the outcome variable directly in the formula i.e. lme(sqrt(y)~x) will return p=1!
#' @param residuals Residuals to be used when constructing the process. Possible values are \code{"individual"} and \code{"cluster"} for \emph{individual} and \emph{cluster-speciffic} residuals, respectively.
#' @param std.type Type of standardization to be used for the residuals when constructing the process.
#' Currently implemeneted options are \code{1} and \code{2} for \eqn{S_i=\hat\sigma^{-1/2}I_{n_i}} and \eqn{S_i=\hat{V}_i^{-1/2}}.
#' @param use.correction.for.imbalance Logical. use \eqn{n_i^{-1/2} S_i} when standardizing the residuals. Defaults to \code{FALSE}.
#' @param subset.fix Two-sided formula. If nonnull, the process \eqn{W^{F^s}} will be constructed using the variables defined on the RHS of the formula. Deafults to \code{NULL} and the process \eqn{W^{F^s}} is not constructed.
#' @param type How to obtain the processes \eqn{W^m}. Possible values are \code{"simulation"} for the simulation approach, \code{"sign.flip"} for the sign-flipping approach and \code{"permutation"} for the permutation approach. When using \code{type="permutation"}, sign-flipping will be used by default if not specified otherwise by the argument \code{force.permutation.with.O}.
#' @param M Number of random simulations/sign-flipps/permutations. Defaults to \code{100}.
#' @param order.by.original Logical. Should the residuals in the the processes \eqn{W^m} be ordered by the original fitted values? Defaults to \code{FALSE}.
#' Makes sense only for \code{type="sign.flip"} and \code{type="permutation"} since when \code{type="simulation"} the ordering is always based on the original predictions.
#' It is programmed such that \eqn{J_i} is reestimated at each iteration \eqn{m}.
#' @param force.permutation.with.O Logical. Should the permutations be used also for the O process? Defaults to \code{FALSE}.
#' @param verbose Logical. Print the current status of the test. Can slow down the algorithm, but it can make it feel faster. Defaults to \code{FALSE}.
#' @return An object of class \code{"gofLMM"} for which \code{plot} and \code{summary} functions are available.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @seealso \code{\link{gof.lmm.pan}}, \code{\link{plot.gofLMM}} and  \code{\link{summary.gofLMM}}
#' @export
#' @examples
#' # simulate some data:
#' N=50
#' set.seed(1)
#' n<-floor(runif(N,min=1,max=15)) #imbalanced
#' betas<-c(1,1,1,15) #don't change! #the last one is only used whe omit.important.predictor=TRUE
#' norm.eps<-FALSE
#' shape=0.5
#' scale=1
#' norm.re.intercept<-FALSE
#' shape.re.intercept=0.5
#' scale.re.intercept=1
#' norm.re.slope<-FALSE
#' shape.re.slope=0.5
#' scale.re.slope=1
#' sim.re.slope=FALSE
#' over.parameterized.model=FALSE #i.e. fit a variable which is not used when generating the data
#' omit.important.predictor=FALSE
#' yy<-NA
#' x22<-NA
#' id<-NA
#' x1<-NA
#' for (gg in 1:N){
#'
#'   id<-c(id,rep(gg,each=n[gg]))
#'   x11<-rep(rbinom(1,size=1,prob=0.4),each=n[gg])
#'   x1<-c(x1,x11)
#'
#'   if (norm.re.intercept==TRUE) re.int<-rnorm(1,sd=sqrt(2)) else re.int<-rgamma(1,shape=shape.re.intercept,scale=scale.re.intercept)-shape.re.intercept*scale.re.intercept
#'
#'   b<-rep(re.int,each=n[gg])
#'
#'   if (norm.re.slope==TRUE) re.slope<-rnorm(1,sd=sqrt(1)) else re.slope<-rgamma(1,shape=shape.re.slope,scale=scale.re.slope)-shape.re.slope*scale.re.slope
#'
#'   b2<-rep(re.slope,each=n[gg])
#'   x2<-1:n[gg]
#'   x4<-runif(n[gg])
#'
#'   if (norm.eps==TRUE) eps<-rnorm(n[gg]) else eps<-rgamma(n[gg],shape=shape,scale=scale)-shape*scale
#'
#'   if (sim.re.slope==TRUE) {
#'     if (omit.important.predictor==FALSE) y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+b2*x2+eps else y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+b2*x2+eps+betas[4]*x4
#'   } else {
#'     if (omit.important.predictor==FALSE) y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+eps else y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+eps+betas[4]*x4
#'   }
#'   yy<-c(yy,y)
#'  x22<-c(x22,x2)
#' }
#' yy<-yy[-1]
#' x22<-x22[-1]
#' x1<-x1[-1]
#' id<-id[-1]
#' x4<-runif(sum(n))
#' aids.art<-data.frame(ptnt=id,outcome=yy,x1=x1,x2=x22,x4=x4)
#' library(nlme)
#' fit<-lme(fixed=outcome~ x2+x1:x2, data=aids.art, random=~x2|ptnt,control=lmeControl( returnObject = TRUE),method="REML" )
#' fit.gof<-gof.lmm(fit,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,subset.fix=outcome~x2,type= "simulation" ,M=25,order.by.original=FALSE,force.permutation.with=FALSE,verbose=TRUE)
#' plot.gofLMM(fit.gof,type=2,subset.M=NULL,xlab="",main="Example")
#' summary.gofLMM(fit.gof)
#'
#' library(nlme)
#' data(Orthodont)
#' Orthodont$Subject<- rep(1:27,each=4)
#' fm1<-lme(distance~age,random=~1|Subject,data=Orthodont,control=lmeControl( returnObject = TRUE),method="REML")
#' gof.fm1<-gof.lmm(fm1,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,subset.fix=NULL,type= "sign.flip" ,M=50,order.by.original=FALSE,force.permutation.with.O=FALSE,verbose=TRUE)
#' plot.gofLMM(gof.fm1,type=2,subset.M=NULL,xlab="",main="Orthodont, model 1")
#' summary.gofLMM(gof.fm1)
#'
#' fm1.1<-lme(distance~age,random=~age|Subject,data=Orthodont,control=lmeControl( returnObject = TRUE),method="REML")
#' gof.fm1.1<-gof.lmm(fm1.1,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,subset.fix=NULL,type= "sign.flip" ,M=50,order.by.original=FALSE,force.permutation.with.O=FALSE,verbose=TRUE)
#' plot.gofLMM(gof.fm1.1 ,type=2,subset.M=NULL,xlab="",main="Orthodont, model 1.1")
#' summary.gofLMM(gof.fm1.1)
#'
#' fm2<-lme(distance~age+Sex,random=~1|Subject,data=Orthodont,control=lmeControl( returnObject = TRUE),method="REML")
#' gof.fm2<-gof.lmm(fm2,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,subset.fix=distance~age,type= "sign.flip" ,M=50,order.by.original=FALSE,force.permutation.with=FALSE,verbose=TRUE)
#' plot.gofLMM(gof.fm2,type=2,subset.M=NULL,xlab="",main="Orthodont, model 2")
#' summary.gofLMM(gof.fm2)
#'
#' fm2.1<-lme(distance~age*Sex,random=~1|Subject,data=Orthodont,control=lmeControl( returnObject = TRUE),method="REML")
#' gof.fm2.1<-gof.lmm(fm2.1,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,subset.fix=NULL,type= "sign.flip" ,M=50,order.by.original=FALSE,force.permutation.with.O=FALSE,verbose=TRUE)
#' plot.gofLMM(gof.fm2.1,type=2,subset.M=NULL,xlab="",main="Orthodont, model 2.1")
#' summary.gofLMM(gof.fm2.1)



gof.lmm<-function(fit,residuals=c("individual","cluster"),std.type=c(1,2),use.correction.for.imbalance=FALSE,subset.fix=NULL,type=c("simulation","sign.flip","permutation"),M=100,order.by.original=FALSE,force.permutation.with.O=FALSE,verbose=FALSE){

####checks, warnings

if (is.null(fit$data)) stop("Model was fitted with keep.data=FALSE. Use keep.data=TRUE.")

if (verbose) cat("Using  \"verbose=FALSE \" slows down the algorithm, but it might feel faster. \n")

if (type=="permutation") cat("type=\"permutation\" is specified. \n Using permutation for the F (and Fs) process, but sign-flipping for O process. \n Get some snack if M is large and model is complex. \n If \"force.permutation.with.O=TRUE\", ignore the warning and so help you god.")


####preliminaries




id<-fit$data[,names(formula(fit$modelStruct$reStr))]

N<-length(unique(id))
n<-table(id)


id.c<-NA
for (ii in 1:N){
id.c<-c(id.c,rep(ii,n[ii]))
}
id.c<-id.c[-1]

if (sum(as.numeric(id)-id.c)!=0) stop("The ID variables needs to be numeric and ordered from 1:N.")

x<-model.matrix(fit, data=fit$data   )

ZZ<- model.matrix(formula(fit$modelStruct$reStr)[[1]],data=fit$data)



###start gof



resI<-residuals(fit, level = 1  )


resP<-residuals(fit, level = 0  )

estI<-fitted(fit,level=1)
 estP<-fitted(fit,level=0)

orI<-order(estI)
orP<-order(estP)



vc<-VarCorr(fit)
sigma.est<-as.numeric(vc[nrow(vc),1])

D<-getVarCov(fit)

beta.f<-fixef(fit)

V<-list()
V.i<-list()
Z<-list()

H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
for (gg in 1:N){
if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
I<-diag(rep(1),n[[gg]])
V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
V.i[[gg]]<-V[[gg]]%^%(-1)
if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
}

H.i<-solve(H)


J<-list()
A<-list()
B<-list()

res.i.c<-resI


for (gg in 1:N){


if (n[gg]!=1) A[[gg]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else A[[gg]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

if (n[gg]!=1) B[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else B[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])



I<-diag(rep(1,n[gg]))

if (residuals=="individual") J[[gg]]<-sigma.est*V.i[[gg]]-(A[[gg]])%*%ginv(B[[gg]])%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]] else J[[gg]]<-I-(A[[gg]]+B[[gg]])%*%ginv(B[[gg]])%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]


if (residuals=="individual") res.i.c[id==gg]<- J[[gg]]%*% resI[id==gg] else  res.i.c[id==gg]<- J[[gg]]%*% resP[id==gg]



}



V.ii.inv<-list()
V.ii<-list()
S.i<-list()

if (residuals=="individual") res.i.c2<-resI else res.i.c2<-resP

respermute<-NA
resIst<-NA
resPst<-NA
for (gg in 1:N){
I<-diag(rep(1,n[gg]))

V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
V.ii[[gg]]<-V[[gg]]%^%(0.5)

resPMp<-matrix(resP[id==gg],ncol=1,nrow=n[gg],byrow=F)
resPMp2<-V.ii.inv[[gg]]%*%resPMp

respermute<-c(respermute,resPMp2)

if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])

resPMpC<-matrix(res.i.c[id==gg],ncol=1,nrow=n[gg],byrow=F)
resPMpC2<-S.i[[gg]]%*%resPMpC
resPMpC2<-resPMpC2

resIst<-c(resIst,resPMpC2)


resPMpCP<-matrix(res.i.c2[id==gg],ncol=1,nrow=n[gg],byrow=F)
resPMpC2P<-S.i[[gg]]%*%resPMpCP
resPMpC2P<-resPMpC2P

resPst<-c(resPst,resPMpC2P)

}

respermute<-respermute[-1]
resIst<-resIst[-1]
 resPst<-resPst[-1]


resoI2<-resIst[orI]
 t01<- estI

for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
ig<-which(round(t01[orI],10)==round(ii,10))
resoI2[ig]<-sum(resoI2[ig])/length(ig)
}

WI2<-1/sqrt(N )*cumsum(resoI2)

resoP2<-resPst[orP]
 t01P<- estP
for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
ig<-which(round(t01P[orP],10)==round(ii,10))
resoP2[ig]<-sum(resoP2[ig])/length(ig)
}

WP2<-1/sqrt(N )*cumsum(resoP2)


##for Fs:
if (!is.null(subset.fix)){
x.subset<-model.matrix(subset.fix, data=fit$data   )
cfs.fix.sub<-fixef(fit)[colnames(x.subset)]

estS<-x.subset%*%cfs.fix.sub
orS<-order(estS)

resoP22<-resPst[orS]
 t01P<- estS
for (ii in as.numeric(names(table(t01P[orS]))[which(table(t01P[orS])>1)])){
ig<-which(round(t01P[orS],10)==round(ii,10))
resoP22[ig]<-sum(resoP22[ig])/length(ig)
}

WP2s<-1/sqrt(N )*cumsum(resoP22)

WsP21<-list()
estSm<-list()
} else {estS<-orS<-WsP21<-estSm<-WP2s<-NULL}




####start sim/sign/permuted proces


if (type=="simulation"){




WsP2<- WsI2 <-list()
estIm<-estPm<-list()

for (iiii in 1:M){

if (verbose) print(paste("Iteration: ",iiii,sep=""))


smp<-rnorm(nrow(x))

newres<-NA
for (gg in 1:N){
newres<-c(newres, V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
}

newres<-newres[-1]



##prvi del procesa

prvi.del.p<-prvi.del<-NA

for (gg in 1:N){

  prvi.del<-c(prvi.del,S.i[[gg]]%*%J[[gg]]%*%(newres[id==gg]))
  if (residuals=="cluster") prvi.del.p<-c(prvi.del.p,S.i[[gg]]%*%(newres[id==gg])) else prvi.del.p<-c(prvi.del.p,sigma.est*S.i[[gg]]%*%V.i[[gg]]%*%(newres[id==gg]))

}

prvi.del<-prvi.del[-1]
prvi.del.p<-prvi.del.p[-1]

prvi.del.o<-prvi.del[orI]

t01<- estI

for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
  ig<-which(round(t01[orI],10)==round(ii,10))
  prvi.del.o[ig]<-sum(prvi.del.o[ig])/length(ig)
}


I<-1/sqrt(N)*cumsum(prvi.del.o)
prvi.del.op<-prvi.del.p[orP]

t01P<- estP
for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
  ig<-which(round(t01P[orP],10)==round(ii,10))
  prvi.del.op[ig]<-sum(prvi.del.op[ig])/length(ig)
}


Ip<-1/sqrt(N)*cumsum(prvi.del.op)


dva.1<-matrix(0,ncol=1,nrow=ncol(x))

for (gg  in 1:N){

  if (n[gg]!=1) dva.1<-dva.1+  t(x[id==gg,])%*%V.i[[gg]]%*%(newres[id==gg]) else dva.1<-dva.1+  matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%(newres[id==gg])

}

drugi.del.p<-drugi.del<-NA

for (gg in 1:N){

  drugi.del<-c(drugi.del,S.i[[gg]]%*%J[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)
  if (residuals=="cluster") drugi.del.p<-c(drugi.del.p,S.i[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1) else drugi.del.p<-c(drugi.del.p,sigma.est*S.i[[gg]]%*%V.i[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)


}

drugi.del<-drugi.del[-1]
drugi.del.p<-drugi.del.p[-1]

drugi.del.o<-drugi.del[orI]


t01<- estI

for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
  ig<-which(round(t01[orI],10)==round(ii,10))
  drugi.del.o[ig]<-sum(drugi.del.o[ig])/length(ig)
}


drugi.del.op<-drugi.del.p[orP]
t01P<- estP
for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
  ig<-which(round(t01P[orP],10)==round(ii,10))
  drugi.del.op[ig]<-sum(drugi.del.op[ig])/length(ig)
}



II<-1/sqrt( N)*cumsum(drugi.del.o)
IIp<-1/sqrt( N)*cumsum(drugi.del.op)

WsI2[[iiii]]<-I-II
WsP2[[iiii]]<-Ip-IIp

estIm[[iiii]]<-estI
estPm[[iiii]]<-estP

if (!is.null(subset.fix)){


##prvi del procesa


prvi.del.op<-prvi.del.p[orS]

t01P<- estS
for (ii in as.numeric(names(table(t01P[orS]))[which(table(t01P[orS])>1)])){
ig<-which(round(t01P[orS],10)==round(ii,10))
prvi.del.op[ig]<-sum(prvi.del.op[ig])/length(ig)
}


Ip<-1/sqrt(N)*cumsum(prvi.del.op)




drugi.del.op<-drugi.del.p[orS]
t01P<- estS
for (ii in as.numeric(names(table(t01P[orS]))[which(table(t01P[orS])>1)])){
ig<-which(round(t01P[orS],10)==round(ii,10))
drugi.del.op[ig]<-sum(drugi.del.op[ig])/length(ig)
}


IIp<-1/sqrt( N)*cumsum(drugi.del.op)

WsP21[[iiii]]<-Ip-IIp
estSm[[iiii]]<-estS
}

}




}

if (type!="simulation"){

if (type=="sign.flip") {

WsP2<- WsI2 <-list()
 estIm<-estPm<-list()

for (iiii in 1:M){

if (verbose) print(paste("Iteration: ",iiii,sep=""))

smp<-sample(c(-1,1),size=sum(n),replace=TRUE)

ys<-NA
for (gg in 1:N){
ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
}
ys<-ys[-1]

datas<-fit$data
datas[,as.character(fit$call$fixed)[2]]<-ys



fits<-suppressWarnings(update(fit,data=datas))



sim.proc<-get.sim.proc(fits, residuals=residuals,std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,subset.fix=subset.fix,order.by.original=order.by.original,or.original.fitted.I=orI,or.original.fitted.P=orP,or.original.fitted.S=orS,
original.fitted.I=estI ,original.fitted.P=estP ,original.fitted.S=estS,
n=n,N=N,x=x,ZZ=ZZ,id=id)

 WsI2[[iiii]]<-sim.proc[[1]]
 WsP2[[iiii]]<-sim.proc[[2]]

if (!is.null(subset.fix)) {

WsP21[[iiii]]<-sim.proc[[3]]
estIm[[iiii]]<-sim.proc[[4]]
estPm[[iiii]]<-sim.proc[[5]]
estSm[[iiii]]<-sim.proc[[6]]
} else {
estIm[[iiii]]<-sim.proc[[3]]
estPm[[iiii]]<-sim.proc[[4]]
}

} #end for

} else { #end if sign.flip

WsP2<- WsI2 <-list()
 estIm<-estPm<-list()

for (iiii in 1:M){

if (verbose) print(paste("Iteration: ",iiii,sep=""))

ys<-NA
for (gg in 1:N){

if (n[gg]==1) smp<-1 else smp<-sample(1:n[gg])
ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute[id==gg])[smp]   )  )
}
ys<-ys[-1]

datas<-fit$data
datas[,as.character(fit$call$fixed)[2]]<-ys



fits<-suppressWarnings(update(fit,data=datas))



sim.procF<-get.sim.proc.F(fits, residuals=residuals,std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,subset.fix=subset.fix,order.by.original=order.by.original,or.original.fitted.P=orP,or.original.fitted.S=orS,
original.fitted.P=estP ,original.fitted.S=estS,
n=n,N=N,x=x,ZZ=ZZ,id=id)


if (!is.null(subset.fix)) {WsP2[[iiii]]<-sim.procF[[1]];WsP21[[iiii]]<-sim.procF[[2]];estPm[[iiii]]<-sim.procF[[3]];estSm[[iiii]]<-sim.procF[[4]]} else  {WsP2[[iiii]]<-sim.procF[[1]];estPm[[iiii]]<-sim.procF[[2]] }


###needed to force sign-flipp for O

if (force.permutation.with.O==FALSE){

smp<-sample(c(-1,1),size=sum(n),replace=TRUE)

ys<-NA
for (gg in 1:N){
ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
}
ys<-ys[-1]

datas<-fit$data
datas[,as.character(fit$call$fixed)[2]]<-ys



fits<-suppressWarnings(update(fit,data=datas))

}

sim.procO<-get.sim.proc.O(fits, residuals=residuals,std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance, order.by.original=order.by.original,or.original.fitted.I=orI,
original.fitted.I=estI ,
n=n,N=N,x=x,ZZ=ZZ,id=id)

WsI2[[iiii]]<-sim.procO[[1]]
estIm[[iiii]]<-sim.procO[[2]]

} #end for

} #end of else

} #end if not sim





res<-list(O=WI2,F=WP2,Om=WsI2,Fm=WsP2,Fs=WP2s,Fsm=WsP21,predO=estI,predOm=estIm,predF=estP,predFm=estPm,predFs=estS,predFsm=estSm)
class(res)<-"gofLMM"
res


} #end of function










#' Plot Function
#'
#' plots the processes which are the result of a call to \code{gof.lmm}
#'
#' @param x an object of class \code{"gofLMM"}, an object returned by a call to \code{\link{gof.lmm}}
#' @param type Type of x-axis. Possible values are 1 for 1:N and 2 for the predicted values. Defaults to 2.
#' @param subset.M How many realizations of \eqn{W^m} should be plotted. Defaults to NULL and all the realizations are plotted.
#' @param ... additional arguments passed to from or to other methods
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @seealso \code{\link{gof.lmm.pan}}, \code{\link{gof.lmm}} and \code{\link{summary.gofLMM}}
#' @export


plot.gofLMM<-function(x,type=2,subset.M=NULL,...){
object<-x
txt1<-expression(W^O)
txt2<-expression(W^F)
txt3<-expression(W^F^S)

if (is.null(object$F)) par(mfrow=c(1,1)) else{ if (is.null(object$Fs)) par(mfrow=c(1,2)) else par(mfrow=c(1,3))}

if (is.null(subset.M)) sbset<-1:length(object$Om) else sbset<-sample(1:length(object$Om),subset.M)


my.plot.gofLMM.part(object$O,object$Om[sbset],type=type,y=object$predO,ym=object$predOm[sbset],ylab=txt1,...)
if (!is.null(object$F)) my.plot.gofLMM.part(object$F,object$Fm[sbset],type=type,y=object$predF,ym=object$predFm[sbset],ylab=txt2,...)

if (!is.null(object$Fs))  my.plot.gofLMM.part(object$Fs,object$Fsm[sbset],type=type,y=object$predFs,ym=object$predFsm[sbset],ylab=txt3,...)



}



#' Summary Function
#'
#' makes a summary of a call to \code{gof.lmm}
#'
#' @param object an object of class \code{"gofLMM"}, an object returned by a call to \code{\link{gof.lmm}}
#' @param ... additional arguments passed to from or to other methods
#' @return a matrix containing KS and CvM test statistics and corresponding \eqn{p}-values for the constructed processes.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @seealso \code{\link{gof.lmm.pan}}, \code{\link{gof.lmm}} and \code{\link{plot.gofLMM}}
#' @export


summary.gofLMM<-function(object,...){

O.s<-test.stat.p.val(object$O,object$Om)
if (!is.null(object$F)) F.s<-test.stat.p.val(object$F,object$Fm) else F.s<-NULL
if (!is.null(object$Fs))  S.s<-test.stat.p.val(object$Fs,object$Fsm) else S.s<-NULL

res<-rbind(O.s,F.s,S.s)
if (is.null(object$F) ) rownames(res)<-paste("O",rownames(res)[1:2],sep=":") else  { if (!is.null(object$Fs)) rownames(res)<-c(paste("O",rownames(res)[1:2],sep=":"),paste("F",rownames(res)[1:2],sep=":")  ,paste("Fs",rownames(res)[1:2],sep=":") ) else rownames(res)<-c(paste("O",rownames(res)[1:2],sep=":"),paste("F",rownames(res)[1:2],sep=":"))}

res

}


#' Print Function
#'
#' prints results from a call to \code{gof.lmm}
#'
#' @param x an object of class \code{"gofLMM"}, an object returned by a call to \code{\link{gof.lmm}}
#' @param ... additional arguments passed to from or to other methods
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @seealso \code{\link{gof.lmm}}, \code{\link{summary.gofLMM}} and \code{\link{plot.gofLMM}}
#' @export


print.gofLMM<-function(x,...){

  cat("Cumsum process.")

}







#' Goodness-of fit test for LMM as proposed by Pan et al.
#'
#' Goodness-of fit test based on cumulative sum stochastic process using the simulation approach proposed by Pan et al.

#'
#' @param fit The result of a call to \code{"nlme"}. The model must be fitted with \code{control=lmeControl( returnObject = TRUE)} and \code{keep.data=TRUE}. An error message is returned otherwise. ID variable must be numeric and ordered from 1:N !
#' @param residuals Residuals to be used when constructing the process. Possible values are \code{"individual"} and \code{"cluster"} for \emph{individual} and \emph{cluster-speciffic} residuals, respectively.
#' @param std.type Type of standardization to be used for the residuals when constructing the process.
#' Currently implemeneted options are \code{1} and \code{2} for \eqn{S_i=\hat\sigma^{-1/2}I_{n_i}} and \eqn{S_i=\hat{V}_i^{-1/2}}.
#' @param use.correction.for.imbalance Logical. use \eqn{n_i^{-1/2} S_i} when standardizing the residuals. Defaults to \code{FALSE}.
#' @param subset.fix Two-sided formula. If nonnull, the process \eqn{W^{F^s}} will be constructed using the variables defined on the RHS of the formula. Deafults to \code{NULL} and the process \eqn{W^{F^s}} is not constructed.
#' @param M Number of random simulations/sign-flipps/permutations. Defaults to \code{100}.
#' @param verbose Logical. Print the current status of the test. Can slow down the algorithm, but it can make it feel faster. Defaults to \code{FALSE}.
#' @return An object of class \code{"gofLMM"} for which \code{plot} and \code{summary} functions are available.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @seealso \code{\link{gof.lmm}}, \code{\link{plot.gofLMM}} and  \code{\link{summary.gofLMM}}
#' @export
#' @examples
#' # simulate some data:
#' N=50
#' set.seed(1)
#' n<-floor(runif(N,min=1,max=15)) #imbalanced
#' betas<-c(1,1,1,15) #don't change! #the last one is only used whe omit.important.predictor=TRUE
#' norm.eps<-FALSE
#' shape=0.5
#' scale=1
#' norm.re.intercept<-FALSE
#' shape.re.intercept=0.5
#' scale.re.intercept=1
#' norm.re.slope<-FALSE
#' shape.re.slope=0.5
#' scale.re.slope=1
#' sim.re.slope=FALSE
#' over.parameterized.model=FALSE #i.e. fit a variable which is not used when generating the data
#' omit.important.predictor=FALSE
#' yy<-NA
#' x22<-NA
#' id<-NA
#' x1<-NA
#' for (gg in 1:N){
#'
#'   id<-c(id,rep(gg,each=n[gg]))
#'   x11<-rep(rbinom(1,size=1,prob=0.4),each=n[gg])
#'   x1<-c(x1,x11)
#'
#'   if (norm.re.intercept==TRUE) re.int<-rnorm(1,sd=sqrt(2)) else re.int<-rgamma(1,shape=shape.re.intercept,scale=scale.re.intercept)-shape.re.intercept*scale.re.intercept
#'
#'   b<-rep(re.int,each=n[gg])
#'
#'   if (norm.re.slope==TRUE) re.slope<-rnorm(1,sd=sqrt(1)) else re.slope<-rgamma(1,shape=shape.re.slope,scale=scale.re.slope)-shape.re.slope*scale.re.slope
#'
#'   b2<-rep(re.slope,each=n[gg])
#'   x2<-1:n[gg]
#'   x4<-runif(n[gg])
#'
#'   if (norm.eps==TRUE) eps<-rnorm(n[gg]) else eps<-rgamma(n[gg],shape=shape,scale=scale)-shape*scale
#'
#'   if (sim.re.slope==TRUE) {
#'     if (omit.important.predictor==FALSE) y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+b2*x2+eps else y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+b2*x2+eps+betas[4]*x4
#'   } else {
#'     if (omit.important.predictor==FALSE) y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+eps else y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+eps+betas[4]*x4
#'   }
#'   yy<-c(yy,y)
#'  x22<-c(x22,x2)
#' }
#' yy<-yy[-1]
#' x22<-x22[-1]
#' x1<-x1[-1]
#' id<-id[-1]
#' x4<-runif(sum(n))
#' aids.art<-data.frame(ptnt=id,outcome=yy,x1=x1,x2=x22,x4=x4)
#' library(nlme)
#' fit<-lme(fixed=outcome~ x2+x1:x2, data=aids.art, random=~x2|ptnt,control=lmeControl( returnObject = TRUE),method="REML" )
#' fit.gof.pan<-gof.lmm.pan(fit,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,subset.fix=outcome~x2,M=25,verbose=TRUE)
#' plot.gofLMM(fit.gof.pan,type=2,subset.M=NULL,xlab="",main="Example")
#' summary.gofLMM(fit.gof.pan)


gof.lmm.pan<-function(fit,residuals=c("individual","cluster"),std.type=c(1,2),use.correction.for.imbalance=FALSE,subset.fix=NULL,M=100,verbose=FALSE){

  ####checks, warnings


  if (verbose) cat("Using  \"verbose=FALSE \" slows down the algorithm, but it might feel faster. \n")


  ####preliminaries
  id<-fit$data[,names(formula(fit$modelStruct$reStr))]


  N<-length(unique(id))
  n<-table(id)



  id.c<-NA
  for (ii in 1:N){
    id.c<-c(id.c,rep(ii,n[ii]))
  }
  id.c<-id.c[-1]

  if (sum(as.numeric(id)-id.c)!=0) stop("The ID variables needs to be numeric and ordered from 1:N.")

  x<-model.matrix(fit, data=fit$data   )

  ZZ<- model.matrix(formula(fit$modelStruct$reStr)[[1]],data=fit$data)



  ###start gof



  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )

  estI<-fitted(fit,level=1)
  estP<-fitted(fit,level=0)

  orI<-order(estI)
  orP<-order(estP)



  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()

  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
  }

  H.i<-solve(H)


  J<-list()
  A<-list()
  B<-list()

  res.i.c<-resI


  for (gg in 1:N){


    if (n[gg]!=1) A[[gg]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else A[[gg]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

    if (n[gg]!=1) B[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else B[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])



    I<-diag(rep(1,n[gg]))

    if (residuals=="individual") J[[gg]]<-sigma.est*V.i[[gg]]-(A[[gg]])%*%ginv(B[[gg]])%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]] else J[[gg]]<-I-(A[[gg]]+B[[gg]])%*%ginv(B[[gg]])%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]


    if (residuals=="individual") res.i.c[id==gg]<- J[[gg]]%*% resI[id==gg] else  res.i.c[id==gg]<- J[[gg]]%*% resP[id==gg]



  }



  V.ii.inv<-list()
  V.ii<-list()
  S.i<-list()
  if (residuals=="individual") res.i.c2<-resI else res.i.c2<-resP

  respermute<-resP
  resIst<-NA
  resPst<-NA
  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
    V.ii[[gg]]<-V[[gg]]%^%(0.5)


    if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])

    resPMpC<-matrix(res.i.c[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-S.i[[gg]]%*%resPMpC
    resPMpC2<-resPMpC2

    resIst<-c(resIst,resPMpC2)


    resPMpCP<-matrix(res.i.c2[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-S.i[[gg]]%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst<-c(resPst,resPMpC2P)

  }

   resIst<-resIst[-1]
  resPst<-resPst[-1]


  resoI2<-resIst[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2<-1/sqrt(N )*cumsum(resoI2)

  resoP2<-resPst[orP]
  t01P<- estP
  for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
    ig<-which(round(t01P[orP],10)==round(ii,10))
    resoP2[ig]<-sum(resoP2[ig])/length(ig)
  }

  WP2<-1/sqrt(N )*cumsum(resoP2)


  ##for Fs:
  if (!is.null(subset.fix)){
    x.subset<-model.matrix(subset.fix, data=fit$data   )
    cfs.fix.sub<-fixef(fit)[colnames(x.subset)]

    estS<-x.subset%*%cfs.fix.sub
    orS<-order(estS)

    resoP22<-resPst[orS]
    t01P<- estS
    for (ii in as.numeric(names(table(t01P[orS]))[which(table(t01P[orS])>1)])){
      ig<-which(round(t01P[orS],10)==round(ii,10))
      resoP22[ig]<-sum(resoP22[ig])/length(ig)
    }

    WP2s<-1/sqrt(N )*cumsum(resoP22)

    WsP21<-list()
    estSm<-list()
  } else {estS<-orS<-WsP21<-estSm<-WP2s<-NULL}




  ####start sim/sign/permuted proces






    WsP2<- WsI2 <-list()
    estIm<-estPm<-list()

    for (iiii in 1:M){

      if (verbose) print(paste("Iteration: ",iiii,sep=""))




      newres<-NA
      for (gg in 1:N){
        smp<-rnorm(1)
        newres<-c(newres,  (respermute*smp)[id==gg])
      }

      newres<-newres[-1]




      ##prvi del procesa

      prvi.del.p<-prvi.del<-NA

      for (gg in 1:N){

        prvi.del<-c(prvi.del,S.i[[gg]]%*%J[[gg]]%*%(newres[id==gg]))
        if (residuals=="cluster") prvi.del.p<-c(prvi.del.p,S.i[[gg]]%*%(newres[id==gg])) else prvi.del.p<-c(prvi.del.p,sigma.est*S.i[[gg]]%*%V.i[[gg]]%*%(newres[id==gg]))

      }

      prvi.del<-prvi.del[-1]
      prvi.del.p<-prvi.del.p[-1]

      prvi.del.o<-prvi.del[orI]

      t01<- estI

      for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
        ig<-which(round(t01[orI],10)==round(ii,10))
        prvi.del.o[ig]<-sum(prvi.del.o[ig])/length(ig)
      }


      I<-1/sqrt(N)*cumsum(prvi.del.o)
      prvi.del.op<-prvi.del.p[orP]

      t01P<- estP
      for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
        ig<-which(round(t01P[orP],10)==round(ii,10))
        prvi.del.op[ig]<-sum(prvi.del.op[ig])/length(ig)
      }


      Ip<-1/sqrt(N)*cumsum(prvi.del.op)


      dva.1<-matrix(0,ncol=1,nrow=ncol(x))

      for (gg  in 1:N){

        if (n[gg]!=1) dva.1<-dva.1+  t(x[id==gg,])%*%V.i[[gg]]%*%(newres[id==gg]) else dva.1<-dva.1+  matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%(newres[id==gg])

      }

      drugi.del.p<-drugi.del<-NA

      for (gg in 1:N){

        drugi.del<-c(drugi.del,S.i[[gg]]%*%J[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)
        if (residuals=="cluster") drugi.del.p<-c(drugi.del.p,S.i[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1) else drugi.del.p<-c(drugi.del.p,sigma.est*S.i[[gg]]%*%V.i[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)


      }

      drugi.del<-drugi.del[-1]
      drugi.del.p<-drugi.del.p[-1]

      drugi.del.o<-drugi.del[orI]


      t01<- estI

      for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
        ig<-which(round(t01[orI],10)==round(ii,10))
        drugi.del.o[ig]<-sum(drugi.del.o[ig])/length(ig)
      }


      drugi.del.op<-drugi.del.p[orP]
      t01P<- estP
      for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
        ig<-which(round(t01P[orP],10)==round(ii,10))
        drugi.del.op[ig]<-sum(drugi.del.op[ig])/length(ig)
      }



      II<-1/sqrt( N)*cumsum(drugi.del.o)
      IIp<-1/sqrt( N)*cumsum(drugi.del.op)

      WsI2[[iiii]]<-I-II
      WsP2[[iiii]]<-Ip-IIp

      estIm[[iiii]]<-estI
      estPm[[iiii]]<-estP

      if (!is.null(subset.fix)){

        ##prvi del procesa


        prvi.del.op<-prvi.del.p[orS]

        t01P<- estS
        for (ii in as.numeric(names(table(t01P[orS]))[which(table(t01P[orS])>1)])){
          ig<-which(round(t01P[orS],10)==round(ii,10))
          prvi.del.op[ig]<-sum(prvi.del.op[ig])/length(ig)
        }


        Ip<-1/sqrt(N)*cumsum(prvi.del.op)




        drugi.del.op<-drugi.del.p[orS]
        t01P<- estS
        for (ii in as.numeric(names(table(t01P[orS]))[which(table(t01P[orS])>1)])){
          ig<-which(round(t01P[orS],10)==round(ii,10))
          drugi.del.op[ig]<-sum(drugi.del.op[ig])/length(ig)
        }


        IIp<-1/sqrt( N)*cumsum(drugi.del.op)

        WsP21[[iiii]]<-Ip-IIp
        estSm[[iiii]]<-estS
      }

    }





  res<-list(O=WI2,F=WP2,Om=WsI2,Fm=WsP2,Fs=WP2s,Fsm=WsP21,predO=estI,predOm=estIm,predF=estP,predFm=estPm,predFs=estS,predFsm=estSm)
  class(res)<-"gofLMM"
  res


} #end of function











#' Internal function
#' @keywords internal


get.sim.proc.fast.ororg<-function(fit, std.type ,use.correction.for.imbalance ,n,N,x,ZZ,id,fittedI, or.fittedI,fittedP,or.fittedP ){


  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )

  estI<- fittedI
  estP<- fittedP

  orI<- or.fittedI
  orP<- or.fittedP


  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()

  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
  }

  H.i<-solve(H)



  res.i.c.ind<-resI

  res.i.c.clust<-resP

  for (gg in 1:N){


    if (n[gg]!=1) A<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else A<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

    if (n[gg]!=1) B<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else B<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])



    I<-diag(rep(1,n[gg]))

    J.ind<-sigma.est*V.i[[gg]]-(A)%*%ginv(B)%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]
    J.clust<-I-(A+B)%*%ginv(B)%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]


    res.i.c.ind[id==gg]<- J.ind%*% resI[id==gg]
    res.i.c.clust[id==gg]<- J.clust%*% resP[id==gg]



  }




  res.i.c2.ind<-resI
  res.i.c2.clust<-resP


  resIst.ind<-NA
  resPst.ind<-NA

  resIst.clust<-NA
  resPst.clust<-NA

  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv<-V[[gg]]%^%(-0.5)




    if (std.type==2) Si<-V.ii.inv else Si<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) Si<-Si/sqrt(n[gg])

    resPMpC<-matrix(res.i.c.ind[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-Si%*%resPMpC
    resPMpC2<-resPMpC2

    resIst.ind<-c(resIst.ind,resPMpC2)


    resPMpCP<-matrix(res.i.c2.ind[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-Si%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst.ind<-c(resPst.ind,resPMpC2P)

    resPMpC<-matrix(res.i.c.clust[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-Si%*%resPMpC
    resPMpC2<-resPMpC2

    resIst.clust<-c(resIst.clust,resPMpC2)


    resPMpCP<-matrix(res.i.c2.clust[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-Si%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst.clust<-c(resPst.clust,resPMpC2P)

  }


  resIst.ind<-resIst.ind[-1]
  resPst.ind<-resPst.ind[-1]

  resIst.clust<-resIst.clust[-1]
  resPst.clust<-resPst.clust[-1]


  resoI2<-resIst.ind[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2.ind<-1/sqrt(N )*cumsum(resoI2)

  resoP2<-resPst.ind[orP]
  t01P<- estP
  for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
    ig<-which(round(t01P[orP],10)==round(ii,10))
    resoP2[ig]<-sum(resoP2[ig])/length(ig)
  }

  WP2.ind<-1/sqrt(N )*cumsum(resoP2)


  resoI2<-resIst.clust[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2.clust<-1/sqrt(N )*cumsum(resoI2)

  resoP2<-resPst.clust[orP]
  t01P<- estP
  for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
    ig<-which(round(t01P[orP],10)==round(ii,10))
    resoP2[ig]<-sum(resoP2[ig])/length(ig)
  }

  WP2.clust<-1/sqrt(N )*cumsum(resoP2)

  list(WI2.ind,WP2.ind,WI2.clust,WP2.clust)

}





#######sim function

#' Goodness-of fit test for LMM, all faster, ordering the residuals by the original fitted values also for the SF and permutation approach
#'
#' Goodness-of fit test based on cumulative sum stochastic process. Used for simulations. Returns only KS and CvM p-values for all 4 methods and individual as well as cluster specific residuals. Fs not implemented here.

#'
#' @param fit The result of a call to \code{"nlme"}. The model must be fitted with \code{control=lmeControl( returnObject = TRUE)} and \code{keep.data=TRUE}; ID variable must be numeric and ordered from 1:N !.
#' @param std.type Type of standardization to be used for the residuals when constructing the process.
#' Currently implemeneted options are \code{1} and \code{2} for \eqn{S_i=\hat\sigma^{-1/2}I_{n_i}} and \eqn{S_i=\hat{V}_i^{-1/2}}.
#' @param use.correction.for.imbalance Logical. use \eqn{n_i^{-1/2} S_i} when standardizing the residuals. Defaults to \code{FALSE}.
#' @param M Number of random simulations/sign-flipps/permutations. Defaults to \code{100}.
#' @param verbose Logical. Print the current status of the test. Can slow down the algorithm, but it can make it feel faster. Defaults to \code{FALSE}.
#' @return KS and CvM pvalues for Pan, Simulation, sign-flip and permutations.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @details for sign.flip and permutation the residuals are ordered by the refitted fitted values.
#' @seealso \code{\link{gof.lmm.pan}} and \code{\link{gof.lmm}}
#' @export
#' @examples
#' # simulate some data:
#' N=50
#' set.seed(1)
#' n<-floor(runif(N,min=1,max=15)) #imbalanced
#' betas<-c(1,1,1,15) #don't change! #the last one is only used whe omit.important.predictor=TRUE
#' norm.eps<-FALSE
#' shape=0.5
#' scale=1
#' norm.re.intercept<-FALSE
#' shape.re.intercept=0.5
#' scale.re.intercept=1
#' norm.re.slope<-FALSE
#' shape.re.slope=0.5
#' scale.re.slope=1
#' sim.re.slope=FALSE
#' over.parameterized.model=FALSE #i.e. fit a variable which is not used when generating the data
#' omit.important.predictor=FALSE
#' yy<-NA
#' x22<-NA
#' id<-NA
#' x1<-NA
#' for (gg in 1:N){
#'
#'   id<-c(id,rep(gg,each=n[gg]))
#'   x11<-rep(rbinom(1,size=1,prob=0.4),each=n[gg])
#'   x1<-c(x1,x11)
#'
#'   if (norm.re.intercept==TRUE) re.int<-rnorm(1,sd=sqrt(2)) else re.int<-rgamma(1,shape=shape.re.intercept,scale=scale.re.intercept)-shape.re.intercept*scale.re.intercept
#'
#'   b<-rep(re.int,each=n[gg])
#'
#'   if (norm.re.slope==TRUE) re.slope<-rnorm(1,sd=sqrt(1)) else re.slope<-rgamma(1,shape=shape.re.slope,scale=scale.re.slope)-shape.re.slope*scale.re.slope
#'
#'   b2<-rep(re.slope,each=n[gg])
#'   x2<-1:n[gg]
#'   x4<-runif(n[gg])
#'
#'   if (norm.eps==TRUE) eps<-rnorm(n[gg]) else eps<-rgamma(n[gg],shape=shape,scale=scale)-shape*scale
#'
#'   if (sim.re.slope==TRUE) {
#'     if (omit.important.predictor==FALSE) y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+b2*x2+eps else y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+b2*x2+eps+betas[4]*x4
#'   } else {
#'     if (omit.important.predictor==FALSE) y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+eps else y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+eps+betas[4]*x4
#'   }
#'   yy<-c(yy,y)
#'  x22<-c(x22,x2)
#' }
#' yy<-yy[-1]
#' x22<-x22[-1]
#' x1<-x1[-1]
#' id<-id[-1]
#' x4<-runif(sum(n))
#' aids.art<-data.frame(ptnt=id,outcome=yy,x1=x1,x2=x22,x4=x4)
#' library(nlme)
#' fit<-lme(fixed=outcome~ x2+x1:x2, data=aids.art, random=~x2|ptnt,control=lmeControl( returnObject = TRUE),method="REML" )
#' gof.lmm.sim.orderbyoriginal(fit,std.type=2,use.correction.for.imbalance=FALSE,M=25,verbose=TRUE)


gof.lmm.sim.orderbyoriginal<-function(fit,std.type=c(1,2),use.correction.for.imbalance=FALSE,M=100,verbose=FALSE){




  id<-fit$data[,names(formula(fit$modelStruct$reStr))]

  N<-length(unique(id))
  n<-table(id)



  x<-model.matrix(fit, data=fit$data   )

  ZZ<- model.matrix(formula(fit$modelStruct$reStr)[[1]],data=fit$data)



  ###start gof



  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )

  estI<-fitted(fit,level=1)
  estP<-fitted(fit,level=0)

  orI<-order(estI)
  orP<-order(estP)



  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()

  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
  }

  H.i<-solve(H)


  J.ind<-list()
  J.clust<-list()
  A<-list()
  B<-list()

  res.i.c.ind<-resI
  res.i.c.clust<-resP

  for (gg in 1:N){


    if (n[gg]!=1) A[[gg]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else A[[gg]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

    if (n[gg]!=1) B[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else B[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])



    I<-diag(rep(1,n[gg]))

    J.ind[[gg]]<-sigma.est*V.i[[gg]]-(A[[gg]])%*%ginv(B[[gg]])%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]
    J.clust[[gg]]<-I-(A[[gg]]+B[[gg]])%*%ginv(B[[gg]])%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]


    res.i.c.ind[id==gg]<- J.ind[[gg]]%*% resI[id==gg]
    res.i.c.clust[id==gg]<- J.clust[[gg]]%*% resP[id==gg]



  }



  V.ii.inv<-list()
  V.ii<-list()
  S.i<-list()

  res.i.c2.ind<-resI
  res.i.c2.clust<-resP

  respermute<-NA
  resIst.ind<-NA
  resPst.ind<-NA
  resIst.clust<-NA
  resPst.clust<-NA
  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
    V.ii[[gg]]<-V[[gg]]%^%(0.5)

    resPMp<-matrix(resP[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMp2<-V.ii.inv[[gg]]%*%resPMp

    respermute<-c(respermute,resPMp2)

    if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])

    resPMpC<-matrix(res.i.c.ind[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-S.i[[gg]]%*%resPMpC
    resPMpC2<-resPMpC2

    resIst.ind<-c(resIst.ind,resPMpC2)

    resPMpC<-matrix(res.i.c.clust[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-S.i[[gg]]%*%resPMpC
    resPMpC2<-resPMpC2

    resIst.clust<-c(resIst.clust,resPMpC2)

    resPMpCP<-matrix(res.i.c2.ind[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-S.i[[gg]]%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst.ind<-c(resPst.ind,resPMpC2P)

    resPMpCP<-matrix(res.i.c2.clust[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-S.i[[gg]]%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst.clust<-c(resPst.clust,resPMpC2P)
  }

  respermute<-respermute[-1]
  resIst.ind<-resIst.ind[-1]
  resPst.ind<-resPst.ind[-1]
  resIst.clust<-resIst.clust[-1]
  resPst.clust<-resPst.clust[-1]

  resoI2<-resIst.ind[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2.ind<-1/sqrt(N )*cumsum(resoI2)

  resoP2<-resPst.ind[orP]
  t01P<- estP
  for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
    ig<-which(round(t01P[orP],10)==round(ii,10))
    resoP2[ig]<-sum(resoP2[ig])/length(ig)
  }

  WP2.ind<-1/sqrt(N )*cumsum(resoP2)



  resoI2<-resIst.clust[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2.clust<-1/sqrt(N )*cumsum(resoI2)

  resoP2<-resPst.clust[orP]
  t01P<- estP
  for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
    ig<-which(round(t01P[orP],10)==round(ii,10))
    resoP2[ig]<-sum(resoP2[ig])/length(ig)
  }

  WP2.clust<-1/sqrt(N )*cumsum(resoP2)



  ####start sim/sign/permuted proces

  ###simulation, Pan approach


  WsP2.ind<- WsI2.ind <-list()

  WsP2.clust<- WsI2.clust <-list()


  for (iiii in 1:M){

    if (verbose) print(paste("Simulation Pan: ",iiii,sep=""))




    newres<-NA
    for (gg in 1:N){
      smp<-rnorm(1)
      newres<-c(newres, ( (resP*smp)[id==gg]))
    }

    newres<-newres[-1]



    ##prvi del procesa

    prvi.del.p.ind<-prvi.del.ind<-prvi.del.p.clust<-prvi.del.clust<-NA

    for (gg in 1:N){

      prvi.del.ind<-c(prvi.del.ind,S.i[[gg]]%*%J.ind[[gg]]%*%(newres[id==gg]))
      prvi.del.clust<-c(prvi.del.clust,S.i[[gg]]%*%J.clust[[gg]]%*%(newres[id==gg]))

      prvi.del.p.clust<-c(prvi.del.p.clust,S.i[[gg]]%*%(newres[id==gg]))
      prvi.del.p.ind<-c(prvi.del.p.ind,sigma.est*S.i[[gg]]%*%V.i[[gg]]%*%(newres[id==gg]))

    }

    prvi.del.ind<-prvi.del.ind[-1]
    prvi.del.p.ind<-prvi.del.p.ind[-1]

    prvi.del.clust<-prvi.del.clust[-1]
    prvi.del.p.clust<-prvi.del.p.clust[-1]

    prvi.del.o<-prvi.del.ind[orI]

    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      prvi.del.o[ig]<-sum(prvi.del.o[ig])/length(ig)
    }


    Iind<-1/sqrt(N)*cumsum(prvi.del.o)

    prvi.del.o<-prvi.del.clust[orI]

    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      prvi.del.o[ig]<-sum(prvi.del.o[ig])/length(ig)
    }


    Iclust<-1/sqrt(N)*cumsum(prvi.del.o)

    prvi.del.op<-prvi.del.p.ind[orP]

    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      prvi.del.op[ig]<-sum(prvi.del.op[ig])/length(ig)
    }


    Ipind<-1/sqrt(N)*cumsum(prvi.del.op)

    prvi.del.op<-prvi.del.p.clust[orP]

    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      prvi.del.op[ig]<-sum(prvi.del.op[ig])/length(ig)
    }


    Ipclust<-1/sqrt(N)*cumsum(prvi.del.op)

    dva.1<-matrix(0,ncol=1,nrow=ncol(x))

    for (gg  in 1:N){

      if (n[gg]!=1) dva.1<-dva.1+  t(x[id==gg,])%*%V.i[[gg]]%*%(newres[id==gg]) else dva.1<-dva.1+  matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%(newres[id==gg])

    }

    drugi.del.p.ind<-drugi.del.ind<-drugi.del.p.clust<-drugi.del.clust<-NA

    for (gg in 1:N){

      drugi.del.ind<-c(drugi.del.ind,S.i[[gg]]%*%J.ind[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)
      drugi.del.clust<-c(drugi.del.clust,S.i[[gg]]%*%J.clust[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)

      drugi.del.p.clust<-c(drugi.del.p.clust,S.i[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)
      drugi.del.p.ind<-c(drugi.del.p.ind,sigma.est*S.i[[gg]]%*%V.i[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)


    }

    drugi.del.ind<-drugi.del.ind[-1]
    drugi.del.p.ind<-drugi.del.p.ind[-1]

    drugi.del.clust<-drugi.del.clust[-1]
    drugi.del.p.clust<-drugi.del.p.clust[-1]

    drugi.del.o<-drugi.del.ind[orI]


    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      drugi.del.o[ig]<-sum(drugi.del.o[ig])/length(ig)
    }


    drugi.del.op<-drugi.del.p.ind[orP]
    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      drugi.del.op[ig]<-sum(drugi.del.op[ig])/length(ig)
    }



    IIind<-1/sqrt( N)*cumsum(drugi.del.o)
    IIpind<-1/sqrt( N)*cumsum(drugi.del.op)

    drugi.del.o<-drugi.del.clust[orI]


    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      drugi.del.o[ig]<-sum(drugi.del.o[ig])/length(ig)
    }


    drugi.del.op<-drugi.del.p.clust[orP]
    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      drugi.del.op[ig]<-sum(drugi.del.op[ig])/length(ig)
    }



    IIclust<-1/sqrt( N)*cumsum(drugi.del.o)
    IIpclust<-1/sqrt( N)*cumsum(drugi.del.op)

    WsI2.ind[[iiii]]<-Iind-IIind
    WsP2.ind[[iiii]]<-Ipind-IIpind

    WsI2.clust[[iiii]]<-Iclust-IIclust
    WsP2.clust[[iiii]]<-Ipclust-IIpclust





  }

  res.sim.ind<-c(
    p.val(  KS(WI2.ind),unlist(lapply(WsI2.ind,KS)) ),
    p.val(  CvM(WI2.ind),unlist(lapply(WsI2.ind,CvM)) ),

    p.val(  KS(WP2.ind),unlist(lapply(WsP2.ind,KS)) ),
    p.val(  CvM(WP2.ind),unlist(lapply(WsP2.ind,CvM)) )
  )

  res.sim.clust<-c(
    p.val(  KS(WI2.clust),unlist(lapply(WsI2.clust,KS)) ),
    p.val(  CvM(WI2.clust),unlist(lapply(WsI2.clust,CvM)) ),

    p.val(  KS(WP2.clust),unlist(lapply(WsP2.clust,KS)) ),
    p.val(  CvM(WP2.clust),unlist(lapply(WsP2.clust,CvM)) )
  )

  ###simulation, our approach




  WsP2.ind<- WsI2.ind <-list()

  WsP2.clust<- WsI2.clust <-list()


  for (iiii in 1:M){

    if (verbose) print(paste("Simulation: ",iiii,sep=""))



    smp<-rnorm(nrow(x))

    newres<-NA
    for (gg in 1:N){
      newres<-c(newres, V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
    }

    newres<-newres[-1]



    ##prvi del procesa

    prvi.del.p.ind<-prvi.del.ind<-prvi.del.p.clust<-prvi.del.clust<-NA

    for (gg in 1:N){

      prvi.del.ind<-c(prvi.del.ind,S.i[[gg]]%*%J.ind[[gg]]%*%(newres[id==gg]))
      prvi.del.clust<-c(prvi.del.clust,S.i[[gg]]%*%J.clust[[gg]]%*%(newres[id==gg]))

      prvi.del.p.clust<-c(prvi.del.p.clust,S.i[[gg]]%*%(newres[id==gg]))
      prvi.del.p.ind<-c(prvi.del.p.ind,sigma.est*S.i[[gg]]%*%V.i[[gg]]%*%(newres[id==gg]))

    }

    prvi.del.ind<-prvi.del.ind[-1]
    prvi.del.p.ind<-prvi.del.p.ind[-1]

    prvi.del.clust<-prvi.del.clust[-1]
    prvi.del.p.clust<-prvi.del.p.clust[-1]

    prvi.del.o<-prvi.del.ind[orI]

    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      prvi.del.o[ig]<-sum(prvi.del.o[ig])/length(ig)
    }


    Iind<-1/sqrt(N)*cumsum(prvi.del.o)

    prvi.del.o<-prvi.del.clust[orI]

    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      prvi.del.o[ig]<-sum(prvi.del.o[ig])/length(ig)
    }


    Iclust<-1/sqrt(N)*cumsum(prvi.del.o)

    prvi.del.op<-prvi.del.p.ind[orP]

    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      prvi.del.op[ig]<-sum(prvi.del.op[ig])/length(ig)
    }


    Ipind<-1/sqrt(N)*cumsum(prvi.del.op)

    prvi.del.op<-prvi.del.p.clust[orP]

    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      prvi.del.op[ig]<-sum(prvi.del.op[ig])/length(ig)
    }


    Ipclust<-1/sqrt(N)*cumsum(prvi.del.op)

    dva.1<-matrix(0,ncol=1,nrow=ncol(x))

    for (gg  in 1:N){

      if (n[gg]!=1) dva.1<-dva.1+  t(x[id==gg,])%*%V.i[[gg]]%*%(newres[id==gg]) else dva.1<-dva.1+  matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%(newres[id==gg])

    }

    drugi.del.p.ind<-drugi.del.ind<-drugi.del.p.clust<-drugi.del.clust<-NA

    for (gg in 1:N){

      drugi.del.ind<-c(drugi.del.ind,S.i[[gg]]%*%J.ind[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)
      drugi.del.clust<-c(drugi.del.clust,S.i[[gg]]%*%J.clust[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)

      drugi.del.p.clust<-c(drugi.del.p.clust,S.i[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)
      drugi.del.p.ind<-c(drugi.del.p.ind,sigma.est*S.i[[gg]]%*%V.i[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)


    }

    drugi.del.ind<-drugi.del.ind[-1]
    drugi.del.p.ind<-drugi.del.p.ind[-1]

    drugi.del.clust<-drugi.del.clust[-1]
    drugi.del.p.clust<-drugi.del.p.clust[-1]

    drugi.del.o<-drugi.del.ind[orI]


    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      drugi.del.o[ig]<-sum(drugi.del.o[ig])/length(ig)
    }


    drugi.del.op<-drugi.del.p.ind[orP]
    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      drugi.del.op[ig]<-sum(drugi.del.op[ig])/length(ig)
    }



    IIind<-1/sqrt( N)*cumsum(drugi.del.o)
    IIpind<-1/sqrt( N)*cumsum(drugi.del.op)

    drugi.del.o<-drugi.del.clust[orI]


    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      drugi.del.o[ig]<-sum(drugi.del.o[ig])/length(ig)
    }


    drugi.del.op<-drugi.del.p.clust[orP]
    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      drugi.del.op[ig]<-sum(drugi.del.op[ig])/length(ig)
    }



    IIclust<-1/sqrt( N)*cumsum(drugi.del.o)
    IIpclust<-1/sqrt( N)*cumsum(drugi.del.op)

    WsI2.ind[[iiii]]<-Iind-IIind
    WsP2.ind[[iiii]]<-Ipind-IIpind

    WsI2.clust[[iiii]]<-Iclust-IIclust
    WsP2.clust[[iiii]]<-Ipclust-IIpclust





  }

  res.sim.our.ind<-c(
    p.val(  KS(WI2.ind),unlist(lapply(WsI2.ind,KS)) ),
    p.val(  CvM(WI2.ind),unlist(lapply(WsI2.ind,CvM)) ),

    p.val(  KS(WP2.ind),unlist(lapply(WsP2.ind,KS)) ),
    p.val(  CvM(WP2.ind),unlist(lapply(WsP2.ind,CvM)) )
  )

  res.sim.our.clust<-c(
    p.val(  KS(WI2.clust),unlist(lapply(WsI2.clust,KS)) ),
    p.val(  CvM(WI2.clust),unlist(lapply(WsI2.clust,CvM)) ),

    p.val(  KS(WP2.clust),unlist(lapply(WsP2.clust,KS)) ),
    p.val(  CvM(WP2.clust),unlist(lapply(WsP2.clust,CvM)) )
  )








  WsP2.ind<- WsI2.ind <-list()
  WsP2.clust<- WsI2.clust <-list()

  for (iiii in 1:M){

    if (verbose) print(paste("Sign-flip: ",iiii,sep=""))

    smp<-sample(c(-1,1),size=sum(n),replace=TRUE)

    ys<-NA
    for (gg in 1:N){
      ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
    }
    ys<-ys[-1]

    datas<-fit$data
    datas[,as.character(fit$call$fixed)[2]]<-ys



    fits<-suppressWarnings(update(fit,data=datas))



    sim.proc<-get.sim.proc.fast.ororg(fits, std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,
                                n=n,N=N,x=x,ZZ=ZZ,id=id,fittedI=estI, or.fittedI=orI,fittedP=estP,or.fittedP=orP)

    WsI2.ind[[iiii]]<-sim.proc[[1]]
    WsP2.ind[[iiii]]<-sim.proc[[2]]


    WsI2.clust[[iiii]]<-sim.proc[[3]]
    WsP2.clust[[iiii]]<-sim.proc[[4]]


  } #end for


  res.sign.ind<-c(
    p.val(  KS(WI2.ind),unlist(lapply(WsI2.ind,KS)) ),
    p.val(  CvM(WI2.ind),unlist(lapply(WsI2.ind,CvM)) ),

    p.val(  KS(WP2.ind),unlist(lapply(WsP2.ind,KS)) ),
    p.val(  CvM(WP2.ind),unlist(lapply(WsP2.ind,CvM)) )
  )
  res.sign.clust<-c(
    p.val(  KS(WI2.clust),unlist(lapply(WsI2.clust,KS)) ),
    p.val(  CvM(WI2.clust),unlist(lapply(WsI2.clust,CvM)) ),

    p.val(  KS(WP2.clust),unlist(lapply(WsP2.clust,KS)) ),
    p.val(  CvM(WP2.clust),unlist(lapply(WsP2.clust,CvM)) )
  )


  WsP2.ind<- WsI2.ind <-list()
  WsP2.clust<- WsI2.clust <-list()


  for (iiii in 1:M){

    if (verbose) print(paste("Permutation: ",iiii,sep=""))

    ys<-NA
    for (gg in 1:N){

      if (n[gg]==1) smp<-1 else smp<-sample(1:n[gg])
      ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute[id==gg])[smp]   )  )
    }
    ys<-ys[-1]

    datas<-fit$data
    datas[,as.character(fit$call$fixed)[2]]<-ys



    fits<-suppressWarnings(update(fit,data=datas))



    sim.proc<-get.sim.proc.fast.ororg(fits, std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,
                                      n=n,N=N,x=x,ZZ=ZZ,id=id,fittedI=estI, or.fittedI=orI,fittedP=estP,or.fittedP=orP)
    WsI2.ind[[iiii]]<-sim.proc[[1]]
    WsP2.ind[[iiii]]<-sim.proc[[2]]


    WsI2.clust[[iiii]]<-sim.proc[[3]]
    WsP2.clust[[iiii]]<-sim.proc[[4]]

  } #end for


  res.perm.ind<-c(
    p.val(  KS(WI2.ind),unlist(lapply(WsI2.ind,KS)) ),
    p.val(  CvM(WI2.ind),unlist(lapply(WsI2.ind,CvM)) ),

    p.val(  KS(WP2.ind),unlist(lapply(WsP2.ind,KS)) ),
    p.val(  CvM(WP2.ind),unlist(lapply(WsP2.ind,CvM)) )
  )
  res.perm.clust<-c(
    p.val(  KS(WI2.clust),unlist(lapply(WsI2.clust,KS)) ),
    p.val(  CvM(WI2.clust),unlist(lapply(WsI2.clust,CvM)) ),

    p.val(  KS(WP2.clust),unlist(lapply(WsP2.clust,KS)) ),
    p.val(  CvM(WP2.clust),unlist(lapply(WsP2.clust,CvM)) )
  )


  res.ind<-c(res.sim.ind,res.sim.our.ind, res.sign.ind,res.perm.ind)
  res.clust<-c(res.sim.clust,res.sim.our.clust, res.sign.clust,res.perm.clust)

  names(res.ind)<-names(res.clust)<-c(
    paste("Sim",c("O.KS","O.CvM","F.KS","F.CvM"),sep=":"),
    paste("SimOur",c("O.KS","O.CvM","F.KS","F.CvM"),sep=":"),
    paste("Sign",c("O.KS","O.CvM","F.KS","F.CvM"),sep=":"),
    paste("Perm",c("O.KS","O.CvM","F.KS","F.CvM"),sep=":")
  )

  resm.ind<-matrix(res.ind,ncol=4,byrow=T)
  resm.clust<-matrix(res.clust,ncol=4,byrow=T)


  colnames(resm.ind)<-colnames(resm.clust)<-c("O.KS","O.CvM","F.KS","F.CvM")
  rownames(resm.ind)<-rownames(resm.clust)<-c("Simulation.Pan","Simulation","sign.flip","permutation")
  list(results.ind=res.ind,results.matrix.ind=resm.ind,results.clust=res.clust,results.matrix.clust=resm.clust)
} #end of function



























#' Internal function
#' @keywords internal


get.sim.proc.fast<-function(fit, std.type ,use.correction.for.imbalance ,n,N,x,ZZ,id ){


  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )

   estI<-fitted(fit,level=1)
   estP<-fitted(fit,level=0)

  orI<-order(estI)
  orP<-order(estP)


  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()

  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
  }

  H.i<-solve(H)



  res.i.c.ind<-resI

  res.i.c.clust<-resP

  for (gg in 1:N){


    if (n[gg]!=1) A<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else A<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

    if (n[gg]!=1) B<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else B<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])



    I<-diag(rep(1,n[gg]))

    J.ind<-sigma.est*V.i[[gg]]-(A)%*%ginv(B)%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]
    J.clust<-I-(A+B)%*%ginv(B)%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]


    res.i.c.ind[id==gg]<- J.ind%*% resI[id==gg]
    res.i.c.clust[id==gg]<- J.clust%*% resP[id==gg]



  }




  res.i.c2.ind<-resI
  res.i.c2.clust<-resP


  resIst.ind<-NA
  resPst.ind<-NA

  resIst.clust<-NA
  resPst.clust<-NA

  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv<-V[[gg]]%^%(-0.5)




    if (std.type==2) Si<-V.ii.inv else Si<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) Si<-Si/sqrt(n[gg])

    resPMpC<-matrix(res.i.c.ind[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-Si%*%resPMpC
    resPMpC2<-resPMpC2

    resIst.ind<-c(resIst.ind,resPMpC2)


    resPMpCP<-matrix(res.i.c2.ind[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-Si%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst.ind<-c(resPst.ind,resPMpC2P)

    resPMpC<-matrix(res.i.c.clust[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-Si%*%resPMpC
    resPMpC2<-resPMpC2

    resIst.clust<-c(resIst.clust,resPMpC2)


    resPMpCP<-matrix(res.i.c2.clust[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-Si%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst.clust<-c(resPst.clust,resPMpC2P)

  }


  resIst.ind<-resIst.ind[-1]
  resPst.ind<-resPst.ind[-1]

  resIst.clust<-resIst.clust[-1]
  resPst.clust<-resPst.clust[-1]


  resoI2<-resIst.ind[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2.ind<-1/sqrt(N )*cumsum(resoI2)

  resoP2<-resPst.ind[orP]
  t01P<- estP
  for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
    ig<-which(round(t01P[orP],10)==round(ii,10))
    resoP2[ig]<-sum(resoP2[ig])/length(ig)
  }

  WP2.ind<-1/sqrt(N )*cumsum(resoP2)


  resoI2<-resIst.clust[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2.clust<-1/sqrt(N )*cumsum(resoI2)

  resoP2<-resPst.clust[orP]
  t01P<- estP
  for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
    ig<-which(round(t01P[orP],10)==round(ii,10))
    resoP2[ig]<-sum(resoP2[ig])/length(ig)
  }

  WP2.clust<-1/sqrt(N )*cumsum(resoP2)

   list(WI2.ind,WP2.ind,WI2.clust,WP2.clust)

}





#######sim function

#' Goodness-of fit test for LMM, all faster
#'
#' Goodness-of fit test based on cumulative sum stochastic process. Used for simulations. Returns only KS and CvM p-values for all 4 methods and individual as well as cluster specific residuals. Fs not implemented here.

#'
#' @param fit The result of a call to \code{"nlme"}. The model must be fitted with \code{control=lmeControl( returnObject = TRUE)} and \code{keep.data=TRUE}; ID variable must be numeric and ordered from 1:N !.
#' @param std.type Type of standardization to be used for the residuals when constructing the process.
#' Currently implemeneted options are \code{1} and \code{2} for \eqn{S_i=\hat\sigma^{-1/2}I_{n_i}} and \eqn{S_i=\hat{V}_i^{-1/2}}.
#' @param use.correction.for.imbalance Logical. use \eqn{n_i^{-1/2} S_i} when standardizing the residuals. Defaults to \code{FALSE}.
#' @param M Number of random simulations/sign-flipps/permutations. Defaults to \code{100}.
#' @param verbose Logical. Print the current status of the test. Can slow down the algorithm, but it can make it feel faster. Defaults to \code{FALSE}.
#' @return KS and CvM pvalues for Pan, Simulation, sign-flip and permutations.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @details for sign.flip and permutation the residuals are ordered by the refitted fitted values.
#' @seealso \code{\link{gof.lmm.pan}} and \code{\link{gof.lmm}}
#' @export
#' @examples
#' # simulate some data:
#' N=50
#' set.seed(1)
#' n<-floor(runif(N,min=1,max=15)) #imbalanced
#' betas<-c(1,1,1,15) #don't change! #the last one is only used whe omit.important.predictor=TRUE
#' norm.eps<-FALSE
#' shape=0.5
#' scale=1
#' norm.re.intercept<-FALSE
#' shape.re.intercept=0.5
#' scale.re.intercept=1
#' norm.re.slope<-FALSE
#' shape.re.slope=0.5
#' scale.re.slope=1
#' sim.re.slope=FALSE
#' over.parameterized.model=FALSE #i.e. fit a variable which is not used when generating the data
#' omit.important.predictor=FALSE
#' yy<-NA
#' x22<-NA
#' id<-NA
#' x1<-NA
#' for (gg in 1:N){
#'
#'   id<-c(id,rep(gg,each=n[gg]))
#'   x11<-rep(rbinom(1,size=1,prob=0.4),each=n[gg])
#'   x1<-c(x1,x11)
#'
#'   if (norm.re.intercept==TRUE) re.int<-rnorm(1,sd=sqrt(2)) else re.int<-rgamma(1,shape=shape.re.intercept,scale=scale.re.intercept)-shape.re.intercept*scale.re.intercept
#'
#'   b<-rep(re.int,each=n[gg])
#'
#'   if (norm.re.slope==TRUE) re.slope<-rnorm(1,sd=sqrt(1)) else re.slope<-rgamma(1,shape=shape.re.slope,scale=scale.re.slope)-shape.re.slope*scale.re.slope
#'
#'   b2<-rep(re.slope,each=n[gg])
#'   x2<-1:n[gg]
#'   x4<-runif(n[gg])
#'
#'   if (norm.eps==TRUE) eps<-rnorm(n[gg]) else eps<-rgamma(n[gg],shape=shape,scale=scale)-shape*scale
#'
#'   if (sim.re.slope==TRUE) {
#'     if (omit.important.predictor==FALSE) y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+b2*x2+eps else y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+b2*x2+eps+betas[4]*x4
#'   } else {
#'     if (omit.important.predictor==FALSE) y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+eps else y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+eps+betas[4]*x4
#'   }
#'   yy<-c(yy,y)
#'  x22<-c(x22,x2)
#' }
#' yy<-yy[-1]
#' x22<-x22[-1]
#' x1<-x1[-1]
#' id<-id[-1]
#' x4<-runif(sum(n))
#' aids.art<-data.frame(ptnt=id,outcome=yy,x1=x1,x2=x22,x4=x4)
#' library(nlme)
#' fit<-lme(fixed=outcome~ x2+x1:x2, data=aids.art, random=~x2|ptnt,control=lmeControl( returnObject = TRUE),method="REML" )
#' gof.lmm.sim(fit,std.type=2,use.correction.for.imbalance=FALSE,M=25,verbose=TRUE)


gof.lmm.sim<-function(fit,std.type=c(1,2),use.correction.for.imbalance=FALSE,M=100,verbose=FALSE){




  id<-fit$data[,names(formula(fit$modelStruct$reStr))]

  N<-length(unique(id))
  n<-table(id)



  x<-model.matrix(fit, data=fit$data   )

  ZZ<- model.matrix(formula(fit$modelStruct$reStr)[[1]],data=fit$data)



  ###start gof



  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )

  estI<-fitted(fit,level=1)
  estP<-fitted(fit,level=0)

  orI<-order(estI)
  orP<-order(estP)



  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()

  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
  }

  H.i<-solve(H)


  J.ind<-list()
  J.clust<-list()
  A<-list()
  B<-list()

  res.i.c.ind<-resI
  res.i.c.clust<-resP

  for (gg in 1:N){


    if (n[gg]!=1) A[[gg]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else A[[gg]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

    if (n[gg]!=1) B[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%t(x[id==gg,]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]]) else B[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- x[id==gg,]%*%H.i%*%matrix(x[id==gg,],ncol=1) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])



    I<-diag(rep(1,n[gg]))

    J.ind[[gg]]<-sigma.est*V.i[[gg]]-(A[[gg]])%*%ginv(B[[gg]])%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]
    J.clust[[gg]]<-I-(A[[gg]]+B[[gg]])%*%ginv(B[[gg]])%*% Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]


    res.i.c.ind[id==gg]<- J.ind[[gg]]%*% resI[id==gg]
    res.i.c.clust[id==gg]<- J.clust[[gg]]%*% resP[id==gg]



  }



  V.ii.inv<-list()
  V.ii<-list()
  S.i<-list()

  res.i.c2.ind<-resI
  res.i.c2.clust<-resP

  respermute<-NA
  resIst.ind<-NA
  resPst.ind<-NA
  resIst.clust<-NA
  resPst.clust<-NA
  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
    V.ii[[gg]]<-V[[gg]]%^%(0.5)

    resPMp<-matrix(resP[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMp2<-V.ii.inv[[gg]]%*%resPMp

    respermute<-c(respermute,resPMp2)

    if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])

    resPMpC<-matrix(res.i.c.ind[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-S.i[[gg]]%*%resPMpC
    resPMpC2<-resPMpC2

    resIst.ind<-c(resIst.ind,resPMpC2)

    resPMpC<-matrix(res.i.c.clust[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-S.i[[gg]]%*%resPMpC
    resPMpC2<-resPMpC2

    resIst.clust<-c(resIst.clust,resPMpC2)

    resPMpCP<-matrix(res.i.c2.ind[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-S.i[[gg]]%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst.ind<-c(resPst.ind,resPMpC2P)

    resPMpCP<-matrix(res.i.c2.clust[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-S.i[[gg]]%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst.clust<-c(resPst.clust,resPMpC2P)
  }

  respermute<-respermute[-1]
  resIst.ind<-resIst.ind[-1]
  resPst.ind<-resPst.ind[-1]
  resIst.clust<-resIst.clust[-1]
  resPst.clust<-resPst.clust[-1]

  resoI2<-resIst.ind[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2.ind<-1/sqrt(N )*cumsum(resoI2)

  resoP2<-resPst.ind[orP]
  t01P<- estP
  for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
    ig<-which(round(t01P[orP],10)==round(ii,10))
    resoP2[ig]<-sum(resoP2[ig])/length(ig)
  }

  WP2.ind<-1/sqrt(N )*cumsum(resoP2)



  resoI2<-resIst.clust[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2.clust<-1/sqrt(N )*cumsum(resoI2)

  resoP2<-resPst.clust[orP]
  t01P<- estP
  for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
    ig<-which(round(t01P[orP],10)==round(ii,10))
    resoP2[ig]<-sum(resoP2[ig])/length(ig)
  }

  WP2.clust<-1/sqrt(N )*cumsum(resoP2)



  ####start sim/sign/permuted proces

  ###simulation, Pan approach


  WsP2.ind<- WsI2.ind <-list()

  WsP2.clust<- WsI2.clust <-list()


  for (iiii in 1:M){

    if (verbose) print(paste("Simulation Pan: ",iiii,sep=""))




    newres<-NA
    for (gg in 1:N){
      smp<-rnorm(1)
      newres<-c(newres, ( (resP*smp)[id==gg]))
    }

    newres<-newres[-1]



    ##prvi del procesa

    prvi.del.p.ind<-prvi.del.ind<-prvi.del.p.clust<-prvi.del.clust<-NA

    for (gg in 1:N){

      prvi.del.ind<-c(prvi.del.ind,S.i[[gg]]%*%J.ind[[gg]]%*%(newres[id==gg]))
      prvi.del.clust<-c(prvi.del.clust,S.i[[gg]]%*%J.clust[[gg]]%*%(newres[id==gg]))

      prvi.del.p.clust<-c(prvi.del.p.clust,S.i[[gg]]%*%(newres[id==gg]))
      prvi.del.p.ind<-c(prvi.del.p.ind,sigma.est*S.i[[gg]]%*%V.i[[gg]]%*%(newres[id==gg]))

    }

    prvi.del.ind<-prvi.del.ind[-1]
    prvi.del.p.ind<-prvi.del.p.ind[-1]

    prvi.del.clust<-prvi.del.clust[-1]
    prvi.del.p.clust<-prvi.del.p.clust[-1]

    prvi.del.o<-prvi.del.ind[orI]

    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      prvi.del.o[ig]<-sum(prvi.del.o[ig])/length(ig)
    }


    Iind<-1/sqrt(N)*cumsum(prvi.del.o)

    prvi.del.o<-prvi.del.clust[orI]

    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      prvi.del.o[ig]<-sum(prvi.del.o[ig])/length(ig)
    }


    Iclust<-1/sqrt(N)*cumsum(prvi.del.o)

    prvi.del.op<-prvi.del.p.ind[orP]

    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      prvi.del.op[ig]<-sum(prvi.del.op[ig])/length(ig)
    }


    Ipind<-1/sqrt(N)*cumsum(prvi.del.op)

    prvi.del.op<-prvi.del.p.clust[orP]

    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      prvi.del.op[ig]<-sum(prvi.del.op[ig])/length(ig)
    }


    Ipclust<-1/sqrt(N)*cumsum(prvi.del.op)

    dva.1<-matrix(0,ncol=1,nrow=ncol(x))

    for (gg  in 1:N){

      if (n[gg]!=1) dva.1<-dva.1+  t(x[id==gg,])%*%V.i[[gg]]%*%(newres[id==gg]) else dva.1<-dva.1+  matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%(newres[id==gg])

    }

    drugi.del.p.ind<-drugi.del.ind<-drugi.del.p.clust<-drugi.del.clust<-NA

    for (gg in 1:N){

      drugi.del.ind<-c(drugi.del.ind,S.i[[gg]]%*%J.ind[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)
      drugi.del.clust<-c(drugi.del.clust,S.i[[gg]]%*%J.clust[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)

      drugi.del.p.clust<-c(drugi.del.p.clust,S.i[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)
      drugi.del.p.ind<-c(drugi.del.p.ind,sigma.est*S.i[[gg]]%*%V.i[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)


    }

    drugi.del.ind<-drugi.del.ind[-1]
    drugi.del.p.ind<-drugi.del.p.ind[-1]

    drugi.del.clust<-drugi.del.clust[-1]
    drugi.del.p.clust<-drugi.del.p.clust[-1]

    drugi.del.o<-drugi.del.ind[orI]


    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      drugi.del.o[ig]<-sum(drugi.del.o[ig])/length(ig)
    }


    drugi.del.op<-drugi.del.p.ind[orP]
    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      drugi.del.op[ig]<-sum(drugi.del.op[ig])/length(ig)
    }



    IIind<-1/sqrt( N)*cumsum(drugi.del.o)
    IIpind<-1/sqrt( N)*cumsum(drugi.del.op)

    drugi.del.o<-drugi.del.clust[orI]


    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      drugi.del.o[ig]<-sum(drugi.del.o[ig])/length(ig)
    }


    drugi.del.op<-drugi.del.p.clust[orP]
    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      drugi.del.op[ig]<-sum(drugi.del.op[ig])/length(ig)
    }



    IIclust<-1/sqrt( N)*cumsum(drugi.del.o)
    IIpclust<-1/sqrt( N)*cumsum(drugi.del.op)

    WsI2.ind[[iiii]]<-Iind-IIind
    WsP2.ind[[iiii]]<-Ipind-IIpind

    WsI2.clust[[iiii]]<-Iclust-IIclust
    WsP2.clust[[iiii]]<-Ipclust-IIpclust





  }

  res.sim.ind<-c(
    p.val(  KS(WI2.ind),unlist(lapply(WsI2.ind,KS)) ),
    p.val(  CvM(WI2.ind),unlist(lapply(WsI2.ind,CvM)) ),

    p.val(  KS(WP2.ind),unlist(lapply(WsP2.ind,KS)) ),
    p.val(  CvM(WP2.ind),unlist(lapply(WsP2.ind,CvM)) )
  )

  res.sim.clust<-c(
    p.val(  KS(WI2.clust),unlist(lapply(WsI2.clust,KS)) ),
    p.val(  CvM(WI2.clust),unlist(lapply(WsI2.clust,CvM)) ),

    p.val(  KS(WP2.clust),unlist(lapply(WsP2.clust,KS)) ),
    p.val(  CvM(WP2.clust),unlist(lapply(WsP2.clust,CvM)) )
  )

###simulation, our approach




  WsP2.ind<- WsI2.ind <-list()

  WsP2.clust<- WsI2.clust <-list()


  for (iiii in 1:M){

    if (verbose) print(paste("Simulation: ",iiii,sep=""))



    smp<-rnorm(nrow(x))

    newres<-NA
    for (gg in 1:N){
      newres<-c(newres, V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
    }

    newres<-newres[-1]



    ##prvi del procesa

    prvi.del.p.ind<-prvi.del.ind<-prvi.del.p.clust<-prvi.del.clust<-NA

    for (gg in 1:N){

      prvi.del.ind<-c(prvi.del.ind,S.i[[gg]]%*%J.ind[[gg]]%*%(newres[id==gg]))
      prvi.del.clust<-c(prvi.del.clust,S.i[[gg]]%*%J.clust[[gg]]%*%(newres[id==gg]))

      prvi.del.p.clust<-c(prvi.del.p.clust,S.i[[gg]]%*%(newres[id==gg]))
      prvi.del.p.ind<-c(prvi.del.p.ind,sigma.est*S.i[[gg]]%*%V.i[[gg]]%*%(newres[id==gg]))

    }

    prvi.del.ind<-prvi.del.ind[-1]
    prvi.del.p.ind<-prvi.del.p.ind[-1]

    prvi.del.clust<-prvi.del.clust[-1]
    prvi.del.p.clust<-prvi.del.p.clust[-1]

    prvi.del.o<-prvi.del.ind[orI]

    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      prvi.del.o[ig]<-sum(prvi.del.o[ig])/length(ig)
    }


    Iind<-1/sqrt(N)*cumsum(prvi.del.o)

    prvi.del.o<-prvi.del.clust[orI]

    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      prvi.del.o[ig]<-sum(prvi.del.o[ig])/length(ig)
    }


    Iclust<-1/sqrt(N)*cumsum(prvi.del.o)

    prvi.del.op<-prvi.del.p.ind[orP]

    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      prvi.del.op[ig]<-sum(prvi.del.op[ig])/length(ig)
    }


    Ipind<-1/sqrt(N)*cumsum(prvi.del.op)

    prvi.del.op<-prvi.del.p.clust[orP]

    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      prvi.del.op[ig]<-sum(prvi.del.op[ig])/length(ig)
    }


    Ipclust<-1/sqrt(N)*cumsum(prvi.del.op)

    dva.1<-matrix(0,ncol=1,nrow=ncol(x))

    for (gg  in 1:N){

      if (n[gg]!=1) dva.1<-dva.1+  t(x[id==gg,])%*%V.i[[gg]]%*%(newres[id==gg]) else dva.1<-dva.1+  matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%(newres[id==gg])

    }

    drugi.del.p.ind<-drugi.del.ind<-drugi.del.p.clust<-drugi.del.clust<-NA

    for (gg in 1:N){

      drugi.del.ind<-c(drugi.del.ind,S.i[[gg]]%*%J.ind[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)
      drugi.del.clust<-c(drugi.del.clust,S.i[[gg]]%*%J.clust[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)

      drugi.del.p.clust<-c(drugi.del.p.clust,S.i[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)
      drugi.del.p.ind<-c(drugi.del.p.ind,sigma.est*S.i[[gg]]%*%V.i[[gg]]%*%x[id==gg,]%*%H.i%*%dva.1)


    }

    drugi.del.ind<-drugi.del.ind[-1]
    drugi.del.p.ind<-drugi.del.p.ind[-1]

    drugi.del.clust<-drugi.del.clust[-1]
    drugi.del.p.clust<-drugi.del.p.clust[-1]

    drugi.del.o<-drugi.del.ind[orI]


    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      drugi.del.o[ig]<-sum(drugi.del.o[ig])/length(ig)
    }


    drugi.del.op<-drugi.del.p.ind[orP]
    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      drugi.del.op[ig]<-sum(drugi.del.op[ig])/length(ig)
    }



    IIind<-1/sqrt( N)*cumsum(drugi.del.o)
    IIpind<-1/sqrt( N)*cumsum(drugi.del.op)

    drugi.del.o<-drugi.del.clust[orI]


    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      drugi.del.o[ig]<-sum(drugi.del.o[ig])/length(ig)
    }


    drugi.del.op<-drugi.del.p.clust[orP]
    t01P<- estP
    for (ii in as.numeric(names(table(t01P[orP]))[which(table(t01P[orP])>1)])){
      ig<-which(round(t01P[orP],10)==round(ii,10))
      drugi.del.op[ig]<-sum(drugi.del.op[ig])/length(ig)
    }



    IIclust<-1/sqrt( N)*cumsum(drugi.del.o)
    IIpclust<-1/sqrt( N)*cumsum(drugi.del.op)

    WsI2.ind[[iiii]]<-Iind-IIind
    WsP2.ind[[iiii]]<-Ipind-IIpind

    WsI2.clust[[iiii]]<-Iclust-IIclust
    WsP2.clust[[iiii]]<-Ipclust-IIpclust





  }

  res.sim.our.ind<-c(
    p.val(  KS(WI2.ind),unlist(lapply(WsI2.ind,KS)) ),
    p.val(  CvM(WI2.ind),unlist(lapply(WsI2.ind,CvM)) ),

    p.val(  KS(WP2.ind),unlist(lapply(WsP2.ind,KS)) ),
    p.val(  CvM(WP2.ind),unlist(lapply(WsP2.ind,CvM)) )
  )

  res.sim.our.clust<-c(
    p.val(  KS(WI2.clust),unlist(lapply(WsI2.clust,KS)) ),
    p.val(  CvM(WI2.clust),unlist(lapply(WsI2.clust,CvM)) ),

    p.val(  KS(WP2.clust),unlist(lapply(WsP2.clust,KS)) ),
    p.val(  CvM(WP2.clust),unlist(lapply(WsP2.clust,CvM)) )
  )








      WsP2.ind<- WsI2.ind <-list()
      WsP2.clust<- WsI2.clust <-list()

      for (iiii in 1:M){

        if (verbose) print(paste("Sign-flip: ",iiii,sep=""))

        smp<-sample(c(-1,1),size=sum(n),replace=TRUE)

        ys<-NA
        for (gg in 1:N){
          ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
        }
        ys<-ys[-1]

        datas<-fit$data
        datas[,as.character(fit$call$fixed)[2]]<-ys



        fits<-suppressWarnings(update(fit,data=datas))



        sim.proc<-get.sim.proc.fast(fits, std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,
                               n=n,N=N,x=x,ZZ=ZZ,id=id)

        WsI2.ind[[iiii]]<-sim.proc[[1]]
        WsP2.ind[[iiii]]<-sim.proc[[2]]


        WsI2.clust[[iiii]]<-sim.proc[[3]]
        WsP2.clust[[iiii]]<-sim.proc[[4]]


      } #end for


      res.sign.ind<-c(
        p.val(  KS(WI2.ind),unlist(lapply(WsI2.ind,KS)) ),
        p.val(  CvM(WI2.ind),unlist(lapply(WsI2.ind,CvM)) ),

        p.val(  KS(WP2.ind),unlist(lapply(WsP2.ind,KS)) ),
        p.val(  CvM(WP2.ind),unlist(lapply(WsP2.ind,CvM)) )
      )
      res.sign.clust<-c(
        p.val(  KS(WI2.clust),unlist(lapply(WsI2.clust,KS)) ),
        p.val(  CvM(WI2.clust),unlist(lapply(WsI2.clust,CvM)) ),

        p.val(  KS(WP2.clust),unlist(lapply(WsP2.clust,KS)) ),
        p.val(  CvM(WP2.clust),unlist(lapply(WsP2.clust,CvM)) )
      )


      WsP2.ind<- WsI2.ind <-list()
      WsP2.clust<- WsI2.clust <-list()


      for (iiii in 1:M){

        if (verbose) print(paste("Permutation: ",iiii,sep=""))

        ys<-NA
        for (gg in 1:N){

          if (n[gg]==1) smp<-1 else smp<-sample(1:n[gg])
          ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute[id==gg])[smp]   )  )
        }
        ys<-ys[-1]

        datas<-fit$data
        datas[,as.character(fit$call$fixed)[2]]<-ys



        fits<-suppressWarnings(update(fit,data=datas))



        sim.proc<-get.sim.proc.fast(fits, std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,
                                    n=n,N=N,x=x,ZZ=ZZ,id=id)

        WsI2.ind[[iiii]]<-sim.proc[[1]]
        WsP2.ind[[iiii]]<-sim.proc[[2]]


        WsI2.clust[[iiii]]<-sim.proc[[3]]
        WsP2.clust[[iiii]]<-sim.proc[[4]]

      } #end for


      res.perm.ind<-c(
        p.val(  KS(WI2.ind),unlist(lapply(WsI2.ind,KS)) ),
        p.val(  CvM(WI2.ind),unlist(lapply(WsI2.ind,CvM)) ),

        p.val(  KS(WP2.ind),unlist(lapply(WsP2.ind,KS)) ),
        p.val(  CvM(WP2.ind),unlist(lapply(WsP2.ind,CvM)) )
      )
      res.perm.clust<-c(
        p.val(  KS(WI2.clust),unlist(lapply(WsI2.clust,KS)) ),
        p.val(  CvM(WI2.clust),unlist(lapply(WsI2.clust,CvM)) ),

        p.val(  KS(WP2.clust),unlist(lapply(WsP2.clust,KS)) ),
        p.val(  CvM(WP2.clust),unlist(lapply(WsP2.clust,CvM)) )
      )


  res.ind<-c(res.sim.ind,res.sim.our.ind, res.sign.ind,res.perm.ind)
  res.clust<-c(res.sim.clust,res.sim.our.clust, res.sign.clust,res.perm.clust)

  names(res.ind)<-names(res.clust)<-c(
  paste("Sim",c("O.KS","O.CvM","F.KS","F.CvM"),sep=":"),
  paste("SimOur",c("O.KS","O.CvM","F.KS","F.CvM"),sep=":"),
  paste("Sign",c("O.KS","O.CvM","F.KS","F.CvM"),sep=":"),
  paste("Perm",c("O.KS","O.CvM","F.KS","F.CvM"),sep=":")
)

 resm.ind<-matrix(res.ind,ncol=4,byrow=T)
 resm.clust<-matrix(res.clust,ncol=4,byrow=T)


 colnames(resm.ind)<-colnames(resm.clust)<-c("O.KS","O.CvM","F.KS","F.CvM")
 rownames(resm.ind)<-rownames(resm.clust)<-c("Simulation.Pan","Simulation","sign.flip","permutation")
 list(results.ind=res.ind,results.matrix.ind=resm.ind,results.clust=res.clust,results.matrix.clust=resm.clust)
} #end of function










#' Internal function
#' @keywords internal


proci<-function(i,res,est){
  or<-order(est[[i]])
  WI<-cumsum( res[[i]][or] )
  t01P<- est[[i]]
  for (ii in as.numeric(names(table(t01P[or]))[which(table(t01P[or])>1)])){
    ig<-which(round(t01P[or],10)==round(ii,10))
    WI[ig]<-sum(WI[ig])/length(ig)
  }
  WI
}

#' Internal function
#' @keywords internal

makeO<-function(res,est,id){

  res.s<-split(res,id)
  est.s<-split(est,id)


  lapply(unique(id),proci,res.s,est.s)


}


#' Internal function
#' @keywords internal


get.sim.proc.i<-function(fit, residuals ,std.type ,use.correction.for.imbalance , order.by.original ,  or.original.fitted.P ,  original.fitted.P ,n,N,x,ZZ,id ){




  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )


  if (order.by.original==TRUE)  estP<-original.fitted.P else  estP<-fitted(fit,level=0)


  if (order.by.original==TRUE)  orP<-or.original.fitted.P  else orP<-order(estP)


  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()

  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
  }

  H.i<-solve(H)




  V.ii.inv<-list()


  if (residuals=="individual") res.i.c2<-resI else res.i.c2<-resP



  resPst<-NA
  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)




    if (std.type==2) Si<-V.ii.inv[[gg]] else Si<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) Si<-Si/sqrt(n[gg])



    resPMpCP<-matrix(res.i.c2[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-Si%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst<-c(resPst,resPMpC2P)

  }



  resPst<-resPst[-1]



  WP2<-makeO(resPst,estP,id)
  WP2



}





#######main function

#' Goodness-of fit test for LMM
#'
#' Goodness-of fit test based on cumulative sum stochastic process for O and John's idea. Not well tested!

#'
#' @param fit The result of a call to \code{"nlme"}. The model must be fitted with \code{control=lmeControl( returnObject = TRUE)} and \code{keep.data=TRUE}. An error message is returned otherwise. ID variable must be numeric and ordered from 1:N ! Canno't use transofrmations of the outcome variable directly in the formula i.e. lme(sqrt(y)~x) will return p=1!
#' @param residuals Residuals to be used when constructing the process. Possible values are \code{"individual"} and \code{"cluster"} for \emph{individual} and \emph{cluster-speciffic} residuals, respectively.
#' @param std.type Type of standardization to be used for the residuals when constructing the process.
#' Currently implemeneted options are \code{1} and \code{2} for \eqn{S_i=\hat\sigma^{-1/2}I_{n_i}} and \eqn{S_i=\hat{V}_i^{-1/2}}.
#' @param use.correction.for.imbalance Logical. use \eqn{n_i^{-1/2} S_i} when standardizing the residuals. Defaults to \code{FALSE}.
#' @param type How to obtain the processes \eqn{W^m}. Possible values are  \code{"sign.flip"} for the sign-flipping approach and \code{"permutation"} for the permutation approach.
#' @param M Number of random simulations/sign-flipps/permutations. Defaults to \code{100}.
#' @param order.by.original Logical. Should the residuals in the the processes \eqn{W^m} be ordered by the original fitted values? Defaults to \code{FALSE}.
#' @param verbose Logical. Print the current status of the test. Can slow down the algorithm, but it can make it feel faster. Defaults to \code{FALSE}.
#' @return An object of class \code{"gofLMM"} for which \code{plot} and \code{summary} functions are available.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @seealso \code{\link{gof.lmm.pan}}, \code{\link{gof.lmm.sim}}
#' @export
#' @examples
#' # simulate some data:
#' N=50
#' set.seed(1)
#' n<-floor(runif(N,min=1,max=15)) #imbalanced
#' betas<-c(1,1,1,15) #don't change! #the last one is only used whe omit.important.predictor=TRUE
#' norm.eps<-FALSE
#' shape=0.5
#' scale=1
#' norm.re.intercept<-FALSE
#' shape.re.intercept=0.5
#' scale.re.intercept=1
#' norm.re.slope<-FALSE
#' shape.re.slope=0.5
#' scale.re.slope=1
#' sim.re.slope=TRUE
#' over.parameterized.model=FALSE #i.e. fit a variable which is not used when generating the data
#' omit.important.predictor=FALSE
#' yy<-NA
#' x22<-NA
#' id<-NA
#' x1<-NA
#' for (gg in 1:N){
#'
#'   id<-c(id,rep(gg,each=n[gg]))
#'   x11<-rep(rbinom(1,size=1,prob=0.4),each=n[gg])
#'   x1<-c(x1,x11)
#'
#'   if (norm.re.intercept==TRUE) re.int<-rnorm(1,sd=sqrt(2)) else re.int<-rgamma(1,shape=shape.re.intercept,scale=scale.re.intercept)-shape.re.intercept*scale.re.intercept
#'
#'   b<-rep(re.int,each=n[gg])
#'
#'   if (norm.re.slope==TRUE) re.slope<-rnorm(1,sd=sqrt(1)) else re.slope<-rgamma(1,shape=shape.re.slope,scale=scale.re.slope)-shape.re.slope*scale.re.slope
#'
#'   b2<-rep(re.slope,each=n[gg])
#'   x2<-1:n[gg]
#'   x4<-runif(n[gg])
#'
#'   if (norm.eps==TRUE) eps<-rnorm(n[gg]) else eps<-rgamma(n[gg],shape=shape,scale=scale)-shape*scale
#'
#'   if (sim.re.slope==TRUE) {
#'     if (omit.important.predictor==FALSE) y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+b2*x2+eps else y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+b2*x2+eps+betas[4]*x4
#'   } else {
#'     if (omit.important.predictor==FALSE) y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+eps else y<-betas[1]+betas[2]*x2+betas[3]*(x11*x2)+b+eps+betas[4]*x4
#'   }
#'   yy<-c(yy,y)
#'  x22<-c(x22,x2)
#' }
#' yy<-yy[-1]
#' x22<-x22[-1]
#' x1<-x1[-1]
#' id<-id[-1]
#' x4<-runif(sum(n))
#' aids.art<-data.frame(ptnt=id,outcome=yy,x1=x1,x2=x22,x4=x4)
#' library(nlme)
#' fit<-lme(fixed=outcome~ x2+x1:x2, data=aids.art, random=~1|ptnt,control=lmeControl( returnObject = TRUE),method="REML" )
#' fit.gof<-gof.lmm.O.type2(fit,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,type= "permutation" ,M=25,order.by.original=FALSE,verbose=TRUE)
#' fit.gof$KS
#' fit2<-lme(fixed=outcome~ x2+x1:x2, data=aids.art, random=~x2|ptnt,control=lmeControl( returnObject = TRUE),method="REML" )
#' fit.gof2<-gof.lmm.O.type2(fit2,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,type= "permutation" ,M=25,order.by.original=FALSE,verbose=TRUE)
#' fit.gof2$KS
#' fit3<-lme(fixed=outcome~ x2+x1:x2, data=aids.art, random=~x1|ptnt,control=lmeControl( returnObject = TRUE),method="REML" )
#' fit.gof3<-gof.lmm.O.type2(fit3,residuals= "individual" ,std.type=2,use.correction.for.imbalance=FALSE,type= "permutation" ,M=25,order.by.original=FALSE,verbose=TRUE)
#' fit.gof3$KS






gof.lmm.O.type2<-function(fit,residuals=c("individual","cluster"),std.type=c(1,2),use.correction.for.imbalance=FALSE,type=c("sign.flip","permutation"),M=100,order.by.original=FALSE,verbose=FALSE){

  ####checks, warnings

  if (is.null(fit$data)) stop("Model was fitted with keep.data=FALSE. Use keep.data=TRUE.")

  if (verbose) cat("Using  \"verbose=FALSE \" slows down the algorithm, but it might feel faster. \n")



  ####preliminaries




  id<-fit$data[,names(formula(fit$modelStruct$reStr))]

  N<-length(unique(id))
  n<-table(id)


  id.c<-NA
  for (ii in 1:N){
    id.c<-c(id.c,rep(ii,n[ii]))
  }
  id.c<-id.c[-1]

  if (sum(as.numeric(id)-id.c)!=0) stop("The ID variables needs to be numeric and ordered from 1:N.")

  x<-model.matrix(fit, data=fit$data   )

  ZZ<- model.matrix(formula(fit$modelStruct$reStr)[[1]],data=fit$data)



  ###start gof



  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )

  estP<-fitted(fit,level=0)

  orP<-order(estP)



  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()

  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
     if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
  }

  H.i<-solve(H)


  if (residuals=="individual") res.i.c2<-resI else res.i.c2<-resP
  V.ii.inv<-list()
  V.ii<-list()
  S.i<-list()

  respermute<-NA
  resPst<-NA
  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
    V.ii[[gg]]<-V[[gg]]%^%(0.5)

    resPMp<-matrix(resP[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMp2<-V.ii.inv[[gg]]%*%resPMp

    respermute<-c(respermute,resPMp2)

    if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])



    resPMpCP<-matrix(res.i.c2[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-S.i[[gg]]%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst<-c(resPst,resPMpC2P)

  }

  respermute<-respermute[-1]
  resPst<-resPst[-1]


  WI2<-makeO(resPst,estP,id)


  ks<-unlist(lapply(WI2,KS))

  cvm<-unlist(lapply(WI2,CvM))



  ksi<-cvmi <-matrix(NA,ncol=length(ks),nrow=M)

WSI<-list()
  for (iiii in 1:M){

    if (verbose) print(paste("Iteration: ",iiii,sep=""))
    if (type=="sign.flip") {
      smp<-sample(c(-1,1),size=sum(n),replace=TRUE)

      ys<-NA
      for (gg in 1:N){
        ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
      }
      ys<-ys[-1]} else {

        ys<-NA
        for (gg in 1:N){

          if (n[gg]==1) smp<-1 else smp<-sample(1:n[gg])
          ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute[id==gg])[smp]   )  )
        }
        ys<-ys[-1]

      }

    datas<-fit$data
    datas[,as.character(fit$call$fixed)[2]]<-ys



    fits<-suppressWarnings(update(fit,data=datas))



    sim.proc<-get.sim.proc.i(fits, residuals=residuals,std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,order.by.original=order.by.original,or.original.fitted.P=orP,original.fitted.P=estP,n=n,N=N,x=x,ZZ=ZZ,id=id)

    ksi[iiii,]<-unlist(lapply(sim.proc,KS))
    cvmi[iiii,]<-unlist(lapply(sim.proc,CvM))
WSI[[iiii]]<-unlist(sim.proc)


  } #end for



  pg.ks<-unlist(lapply(1:N,function(i,x,y) p.val(x[[i]],y[,i]),ks,ksi  ))
  pg.cvm<-unlist(lapply(1:N,function(i,x,y) p.val(x[[i]],y[,i]),cvm,cvmi  ))

  ts.ks<--2*sum(log(pg.ks))
  ts.cvm<--2*sum(log(pg.cvm))

  p.ks<-pchisq( ts.ks, df=  2*N,lower.tail=F )
  p.cvm<-pchisq( ts.cvm, df=  2*N,lower.tail=F )


  res<-list(KS=p.ks,CvM=p.cvm,WI=unlist(WI2),WIm=WSI)
  res


} #end of function





#' Internal function
#' @keywords internal




gof.lmm.O.type2.i<-function(fit,residuals ,std.type,use.correction.for.imbalance,type,M,order.by.original,id,N,n,ZZ,x){



  ###start gof



  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )

  estP<-fitted(fit,level=0)

  orP<-order(estP)



  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()

  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
  }

  H.i<-solve(H)


  if (residuals=="individual") res.i.c2<-resI else res.i.c2<-resP
  V.ii.inv<-list()
  V.ii<-list()
  S.i<-list()

  respermute<-NA
  resPst<-NA
  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
    V.ii[[gg]]<-V[[gg]]%^%(0.5)

    resPMp<-matrix(resP[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMp2<-V.ii.inv[[gg]]%*%resPMp

    respermute<-c(respermute,resPMp2)

    if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])



    resPMpCP<-matrix(res.i.c2[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-S.i[[gg]]%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst<-c(resPst,resPMpC2P)

  }

  respermute<-respermute[-1]
  resPst<-resPst[-1]


  WI2<-makeO(resPst,estP,id)


  ks<-unlist(lapply(WI2,KS))

  cvm<-unlist(lapply(WI2,CvM))



  ksi<-cvmi <-matrix(NA,ncol=length(ks),nrow=M)

  WSI<-list()
  for (iiii in 1:M){

     if (type=="sign.flip") {
      smp<-sample(c(-1,1),size=sum(n),replace=TRUE)

      ys<-NA
      for (gg in 1:N){
        ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
      }
      ys<-ys[-1]} else {

        ys<-NA
        for (gg in 1:N){

          if (n[gg]==1) smp<-1 else smp<-sample(1:n[gg])
          ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute[id==gg])[smp]   )  )
        }
        ys<-ys[-1]

      }

    datas<-fit$data
    datas[,as.character(fit$call$fixed)[2]]<-ys



    fits<-suppressWarnings(update(fit,data=datas))



    sim.proc<-get.sim.proc.i(fits, residuals=residuals,std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,order.by.original=order.by.original,or.original.fitted.P=orP,original.fitted.P=estP,n=n,N=N,x=x,ZZ=ZZ,id=id)

    ksi[iiii,]<-unlist(lapply(sim.proc,KS))
    cvmi[iiii,]<-unlist(lapply(sim.proc,CvM))
    WSI[[iiii]]<-unlist(sim.proc)


  } #end for



  pg.ks<-unlist(lapply(1:N,function(i,x,y) p.val(x[[i]],y[,i]),ks,ksi  ))
  pg.cvm<-unlist(lapply(1:N,function(i,x,y) p.val(x[[i]],y[,i]),cvm,cvmi  ))

  ts.ks<--2*sum(log(pg.ks))
  ts.cvm<--2*sum(log(pg.cvm))

  c(ts.ks,ts.cvm)


} #end of function



#' Goodness-of fit test for LMM
#'
#' Goodness-of fit test based on cumulative sum stochastic process for Oi and John's idea using full parametric bootstrap to obtain p-values. Not well tested!

#'
#' @param fit The result of a call to \code{"nlme"}. The model must be fitted with \code{control=lmeControl( returnObject = TRUE)} and \code{keep.data=TRUE}. An error message is returned otherwise. ID variable must be numeric and ordered from 1:N ! Canno't use transofrmations of the outcome variable directly in the formula i.e. lme(sqrt(y)~x) will return p=1!
#' @param residuals Residuals to be used when constructing the process. Possible values are \code{"individual"} and \code{"cluster"} for \emph{individual} and \emph{cluster-speciffic} residuals, respectively.
#' @param std.type Type of standardization to be used for the residuals when constructing the process.
#' Currently implemeneted options are \code{1} and \code{2} for \eqn{S_i=\hat\sigma^{-1/2}I_{n_i}} and \eqn{S_i=\hat{V}_i^{-1/2}}.
#' @param use.correction.for.imbalance Logical. use \eqn{n_i^{-1/2} S_i} when standardizing the residuals. Defaults to \code{FALSE}.
#' @param type How to obtain the processes \eqn{W^m}. Possible values are  \code{"sign.flip"} for the sign-flipping approach and \code{"permutation"} for the permutation approach.
#' @param M Number of random simulations/sign-flipps/permutations. Defaults to \code{100}.
#' @param B Number of boot replications. Defaults to \code{100}.
#' @param order.by.original Logical. Should the residuals in the the processes \eqn{W^m} be ordered by the original fitted values? Defaults to \code{FALSE}.
#' @param verbose Logical. Print the current status of the test. Can slow down the algorithm, but it can make it feel faster. Defaults to \code{FALSE}.
#' @return An object of class \code{"gofLMM"} for which \code{plot} and \code{summary} functions are available.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @seealso \code{\link{gof.lmm.pan}}, \code{\link{gof.lmm.sim}}
#' @export


gof.lmm.O.type2.boot<-function(fit,residuals=c("individual","cluster"),std.type=c(1,2),use.correction.for.imbalance=FALSE,type=c("sign.flip","permutation"),M=100,B=100,order.by.original=FALSE,verbose=FALSE){

  ####checks, warnings

  if (is.null(fit$data)) stop("Model was fitted with keep.data=FALSE. Use keep.data=TRUE.")

  if (verbose) cat("Using  \"verbose=FALSE \" slows down the algorithm, but it might feel faster. Get some snack as this might take a while. \n")



  ####preliminaries




  id<-fit$data[,names(formula(fit$modelStruct$reStr))]

  N<-length(unique(id))
  n<-table(id)


  id.c<-NA
  for (ii in 1:N){
    id.c<-c(id.c,rep(ii,n[ii]))
  }
  id.c<-id.c[-1]

  if (sum(as.numeric(id)-id.c)!=0) stop("The ID variables needs to be numeric and ordered from 1:N.")

  x<-model.matrix(fit, data=fit$data   )

  ZZ<- model.matrix(formula(fit$modelStruct$reStr)[[1]],data=fit$data)



  ###start gof



  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )

  estP<-fitted(fit,level=0)

  orP<-order(estP)



  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()

  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
  }

  H.i<-solve(H)


  if (residuals=="individual") res.i.c2<-resI else res.i.c2<-resP
  V.ii.inv<-list()
  V.ii<-list()
  S.i<-list()

  respermute<-NA
  resPst<-NA
  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
    V.ii[[gg]]<-V[[gg]]%^%(0.5)

    resPMp<-matrix(resP[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMp2<-V.ii.inv[[gg]]%*%resPMp

    respermute<-c(respermute,resPMp2)

    if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])



    resPMpCP<-matrix(res.i.c2[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2P<-S.i[[gg]]%*%resPMpCP
    resPMpC2P<-resPMpC2P

    resPst<-c(resPst,resPMpC2P)

  }

  respermute<-respermute[-1]
  resPst<-resPst[-1]


  WI2<-makeO(resPst,estP,id)


  ks<-unlist(lapply(WI2,KS))

  cvm<-unlist(lapply(WI2,CvM))



  ksi<-cvmi <-matrix(NA,ncol=length(ks),nrow=M)

  WSI<-list()
  for (iiii in 1:M){

     if (type=="sign.flip") {
      smp<-sample(c(-1,1),size=sum(n),replace=TRUE)

      ys<-NA
      for (gg in 1:N){
        ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
      }
      ys<-ys[-1]} else {

        ys<-NA
        for (gg in 1:N){

          if (n[gg]==1) smp<-1 else smp<-sample(1:n[gg])
          ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute[id==gg])[smp]   )  )
        }
        ys<-ys[-1]

      }

    datas<-fit$data
    datas[,as.character(fit$call$fixed)[2]]<-ys



    fits<-suppressWarnings(update(fit,data=datas))



    sim.proc<-get.sim.proc.i(fits, residuals=residuals,std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,order.by.original=order.by.original,or.original.fitted.P=orP,original.fitted.P=estP,n=n,N=N,x=x,ZZ=ZZ,id=id)

    ksi[iiii,]<-unlist(lapply(sim.proc,KS))
    cvmi[iiii,]<-unlist(lapply(sim.proc,CvM))
    WSI[[iiii]]<-unlist(sim.proc)


  } #end for



  pg.ks<-unlist(lapply(1:N,function(i,x,y) p.val(x[[i]],y[,i]),ks,ksi  ))
  pg.cvm<-unlist(lapply(1:N,function(i,x,y) p.val(x[[i]],y[,i]),cvm,cvmi  ))

  ts.ks<--2*sum(log(pg.ks))
  ts.cvm<--2*sum(log(pg.cvm))

  res.boot<-matrix(NA,ncol=2,nrow=B)

  for (ii in 1:B){
    if (verbose) print(paste("Bootstrap Iteration: ",ii,sep=""))

ys<-NA
    for (jj in 1:N){
     mui<-x[id==jj,]%*%matrix(fixef(fit),ncol=1)
     vari<-Z[[jj]]%*%D%*%t(Z[[jj]])+sigma.est*diag(rep(1,n[jj]))
     ys<-c(ys,rmvnorm(1,mui,vari))
    }
ys<-ys[-1]
datas<-fit$data
datas[,as.character(fit$call$fixed)[2]]<-ys



fits<-suppressWarnings(update(fit,data=datas))
res.boot[ii,]<-gof.lmm.O.type2.i(fit=fits,residuals=residuals ,std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,type=type,M=M,
                  order.by.original=order.by.original,id=id,N=N,n=n,ZZ=ZZ,x=x)


  }


p.ks<-p.val(ts.ks,res.boot[,1])
p.cvm<-p.val(ts.cvm,res.boot[,2])

list(p.ks=p.ks,p.cvm=p.cvm,ts.ks=ts.ks,ts.cvm=ts.cvm)

} #end of function








#' Internal function
#' @keywords internal

get.sim.proc.O.test<-function(fit, residuals ,std.type ,use.correction.for.imbalance ,order.by.original,n,N,x,ZZ,id, est.original,or.original ){



resI<-residuals(fit, level = 1  )


resP<-residuals(fit, level = 0  )

estI<-fitted(fit,level=1)
estP<-fitted(fit,level=0)

orI<-order(estI)



vc<-VarCorr(fit)
sigma.est<-as.numeric(vc[nrow(vc),1])

D<-getVarCov(fit)

beta.f<-fixef(fit)

V<-list()
V.i<-list()
Z<-list()
Xi<-list()
Zb<-list()
H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
for (gg in 1:N){
if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
I<-diag(rep(1),n[[gg]])
V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
V.i[[gg]]<-V[[gg]]%^%(-1)
if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
if (n[gg]!=1) Xi[[gg]]<-x[id==gg,] else Xi[[gg]]<-matrix(x[id==gg,],nrow=1)

}

H.i<-solve(H)


#A<-list()
#B<-list()

res.i.c<-resI

#mm=0
#for (gg in 1:N){
#  for (jj in 1:N){
#    mm=mm+1


#    if (jj==gg){

#      zdz<-  Z[[gg]]%*%D%*%t(Z[[gg]])
#      cpd<-Xi[[gg]]%*%H.i%*%t(Xi[[gg]])

      ###A[[mm]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- Xi[[gg]]%*%H.i%*%t(Xi[[gg]]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])
      ###B[[mm]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- Xi[[gg]]%*%H.i%*%t(Xi[[gg]]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

#      vzd<-V.i[[gg]]%*%( V[[gg]]- cpd )%*%V.i[[gg]]%*%zdz

#      A[[mm]]<-sigma.est*vzd
#      B[[mm]]<-zdz %*%vzd
#    } else {

#      zdzj<-Z[[jj]]%*%D%*%t(Z[[jj]])
#      zdzg<-Z[[gg]]%*%D%*%t(Z[[gg]])
#      cpdj<-Xi[[gg]]%*%H.i%*%t(Xi[[jj]])

      #####A[[mm]]<--sigma.est*V.i[[gg]]%*%(  Xi[[gg]]%*%H.i%*%t(Xi[[jj]]) )%*%V.i[[jj]]%*%Z[[jj]]%*%D%*%t(Z[[jj]])
      ####B[[mm]]<--Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%(  Xi[[gg]]%*%H.i%*%t(Xi[[jj]]) )%*%V.i[[jj]]%*%Z[[jj]]%*%D%*%t(Z[[jj]])

#      vzdj<-V.i[[gg]]%*%(  cpdj )%*%V.i[[jj]]%*%zdzj

#      A[[mm]]<--sigma.est*vzdj
#      B[[mm]]<--zdzg %*%vzdj


#    }

#  }}

#aa<-bb<-list()
#mm=0

#for (gg in 1:N){

#for (jj in 1:N){
#mm=mm+1
#if (jj==1) aa[[gg]]<-A[[mm]] else aa[[gg]]<-cbind(aa[[gg]],A[[mm]])
#if (jj==1) bb[[gg]]<-B[[mm]] else bb[[gg]]<-cbind(bb[[gg]],B[[mm]])
#}}

#for (gg in 1:N){
#if (gg==1) AA<-aa[[gg]] else AA<-rbind(AA,aa[[gg]])
#if (gg==1) BB<-bb[[gg]] else BB<-rbind(BB,bb[[gg]])

#}


NN<-nrow(x)

Vd<-matrix(0,NN,NN)
for (i in 1:N){

  is<-which(id==i)
  Vd[ is,is]<-V[[i]]

}

IN<-diag(rep(1,NN))
Vdi<-solve(Vd)
sVdi<-sigma.est*Vdi
p1<-(IN-sVdi)
p2<-( Vd-x%*%H.i%*%t(x)   )%*%p1
AA<-sVdi%*%p2
BB<-p1%*%p2
Zb<-matrix(estI-estP,ncol=1)


#if (residuals=="individual") res.i.c<-resI-AA%*%ginv(BB)%*%Zb else res.i.c<-resP-(AA+BB)%*%ginv(BB)%*%Zb
if (residuals=="individual") res.i.c<-resI-AA%*%my.MP(BB)%*%Zb else res.i.c<-resP-(AA+BB)%*%my.MP(BB)%*%Zb




V.ii.inv<-list()
V.ii<-list()
S.i<-list()




resIst<-NA


for (gg in 1:N){
I<-diag(rep(1,n[gg]))

V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
V.ii[[gg]]<-V[[gg]]%^%(0.5)



if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])

resPMpC<-matrix(res.i.c[id==gg],ncol=1,nrow=n[gg],byrow=F)
resPMpC2<-S.i[[gg]]%*%resPMpC
resPMpC2<-resPMpC2

resIst<-c(resIst,resPMpC2)


}



resIst<-resIst[-1]


if (order.by.original==TRUE) {estI=est.original  ; orI= or.original }

resoI2<-resIst[orI]
 t01<- estI

for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
ig<-which(round(t01[orI],10)==round(ii,10))
resoI2[ig]<-sum(resoI2[ig])/length(ig)
}

WI2<-1/sqrt(N )*cumsum(resoI2)

list(WI2,estI)

}


#' Internal function
#' @keywords internal


my.MP <- function(A, eps=10^(-8))   {

 PV<-eigen(A,symmetric=T)
 V0<-IV<-PV$values
 IV[abs(V0)>eps]<-1/V0[abs(V0)>eps]
 IV[abs(V0)<=eps]<-0
 Ainv<-PV$vectors%*%(IV*t(PV$vectors  )  )
return(Ainv)


}




#' Goodness-of fit test for LMM
#'
#' Goodness-of fit test based on cumulative sum stochastic process for O using non-diagonal blocked matrices A and B. An error occurs often when calculating the MP generalized inverse of the matrix B, which is due to Lapack routine. Can be very slow and inefficient when n and ni are large. Now I replaced the ginv() from MASS by my.MP which is the MP inverse as suggested by Demidenko p.51 but the error persits, so it must occur in the fitting of lme.

#'
#' @param fit The result of a call to \code{"nlme"}. The model must be fitted with \code{control=lmeControl( returnObject = TRUE)} and \code{keep.data=TRUE}. An error message is returned otherwise. ID variable must be numeric and ordered from 1:N ! Cannot use transformations of the outcome variable directly in the formula i.e. lme(sqrt(y)~x) will return p=1!
#' @param residuals Residuals to be used when constructing the process.
#' @param std.type Type of standardization to be used for the residuals when constructing the process.
#' Currently implemeneted options are \code{1} and \code{2} for \eqn{S_i=\hat\sigma^{-1/2}I_{n_i}} and \eqn{S_i=\hat{V}_i^{-1/2}}.
#' @param use.correction.for.imbalance Logical. use \eqn{n_i^{-1/2} S_i} when standardizing the residuals. Defaults to \code{FALSE}.
#' @param type How to obtain the processes \eqn{W^m}. Possible values are  \code{"sign.flip"} for the sign-flipping approach and \code{"permutation"} for the permutation approach.
#' @param M Number of random simulations/sign-flipps/permutations. Defaults to \code{100}.
#' @param order.by.original Order the residuals by original fitted values? Defaults to FALSE.
#' @param verbose Logical. Print the current status of the test. Can slow down the algorithm, but it can make it feel faster. Defaults to \code{FALSE}.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @seealso \code{\link{gof.lmm}}
#' @export


gof.lmm.O.test<-function(fit,residuals="individual",std.type=c(1,2),use.correction.for.imbalance=FALSE,type=c("sign.flip","permutation"),M=100,order.by.original=FALSE,verbose=FALSE){

####checks, warnings

if (is.null(fit$data)) stop("Model was fitted with keep.data=FALSE. Use keep.data=TRUE.")

if (verbose) cat("Using  \"verbose=TRUE \" slows down the algorithm, but it might feel faster. \n")



####preliminaries




id<-fit$data[,names(formula(fit$modelStruct$reStr))]

N<-length(unique(id))
n<-table(id)


id.c<-NA
for (ii in 1:N){
id.c<-c(id.c,rep(ii,n[ii]))
}
id.c<-id.c[-1]

if (sum(as.numeric(id)-id.c)!=0) stop("The ID variables needs to be numeric and ordered from 1:N.")

x<-model.matrix(fit, data=fit$data   )

ZZ<- model.matrix(formula(fit$modelStruct$reStr)[[1]],data=fit$data)



###start gof



resI<-residuals(fit, level = 1  )


resP<-residuals(fit, level = 0  )

estI<-fitted(fit,level=1)
 estP<-fitted(fit,level=0)

orI<-order(estI)
orP<-order(estP)



vc<-VarCorr(fit)
sigma.est<-as.numeric(vc[nrow(vc),1])

D<-getVarCov(fit)

beta.f<-fixef(fit)

V<-list()
V.i<-list()
Z<-list()
Xi<-list()
Zb<-list()
H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
for (gg in 1:N){
if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
I<-diag(rep(1),n[[gg]])
V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
V.i[[gg]]<-V[[gg]]%^%(-1)
if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
if (n[gg]!=1) Xi[[gg]]<-x[id==gg,] else Xi[[gg]]<-matrix(x[id==gg,],nrow=1)

}

H.i<-solve(H)


#A<-list()
#B<-list()

res.i.c<-resI

#mm=0
#for (gg in 1:N){
#for (jj in 1:N){
#mm=mm+1


#if (jj==gg){

 # zdz<-  Z[[gg]]%*%D%*%t(Z[[gg]])
#  cpd<-Xi[[gg]]%*%H.i%*%t(Xi[[gg]])

      ###A[[mm]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- Xi[[gg]]%*%H.i%*%t(Xi[[gg]]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])
      ####B[[mm]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- Xi[[gg]]%*%H.i%*%t(Xi[[gg]]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

 # vzd<-V.i[[gg]]%*%( V[[gg]]- cpd )%*%V.i[[gg]]%*%zdz

  #A[[mm]]<-sigma.est*vzd
  #B[[mm]]<-zdz %*%vzd
  #} else {

   # zdzj<-Z[[jj]]%*%D%*%t(Z[[jj]])
  #  zdzg<-Z[[gg]]%*%D%*%t(Z[[gg]])
  #  cpdj<-Xi[[gg]]%*%H.i%*%t(Xi[[jj]])

      ####A[[mm]]<--sigma.est*V.i[[gg]]%*%(  Xi[[gg]]%*%H.i%*%t(Xi[[jj]]) )%*%V.i[[jj]]%*%Z[[jj]]%*%D%*%t(Z[[jj]])
      ####B[[mm]]<--Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%(  Xi[[gg]]%*%H.i%*%t(Xi[[jj]]) )%*%V.i[[jj]]%*%Z[[jj]]%*%D%*%t(Z[[jj]])

#vzdj<-V.i[[gg]]%*%(  cpdj )%*%V.i[[jj]]%*%zdzj

#A[[mm]]<--sigma.est*vzdj
#B[[mm]]<--zdzg %*%vzdj


#}

#}}

#aa<-bb<-list()
#mm=0

#for (gg in 1:N){

#for (jj in 1:N){
#mm=mm+1
#if (jj==1) aa[[gg]]<-A[[mm]] else aa[[gg]]<-cbind(aa[[gg]],A[[mm]])
#if (jj==1) bb[[gg]]<-B[[mm]] else bb[[gg]]<-cbind(bb[[gg]],B[[mm]])
#}}

#for (gg in 1:N){
#if (gg==1) AA<-aa[[gg]] else AA<-rbind(AA,aa[[gg]])
#if (gg==1) BB<-bb[[gg]] else BB<-rbind(BB,bb[[gg]])

#}

NN<-nrow(x)

Vd<-Vdi<-matrix(0,NN,NN)
for (i in 1:N){

  is<-which(id==i)
  Vd[ is,is]<-V[[i]]
  Vdi[ is,is]<-V.i[[i]]
}

IN<-diag(rep(1,NN))
#Vdi<-solve(Vd)
sVdi<-sigma.est*Vdi
p1<-(IN-sVdi)
p2<-( Vd-x%*%H.i%*%t(x)   )%*%p1
AA<-sVdi%*%p2
BB<-p1%*%p2





Zb<-matrix(estI-estP,ncol=1)

#if (residuals=="individual") res.i.c<-resI-AA%*%ginv(BB)%*%Zb else res.i.c<-resP-(AA+BB)%*%ginv(BB)%*%Zb
if (residuals=="individual") res.i.c<-resI-AA%*%my.MP(BB)%*%Zb else res.i.c<-resP-(AA+BB)%*%my.MP(BB)%*%Zb



V.ii.inv<-list()
V.ii<-list()
S.i<-list()


respermute<-NA
resIst<-NA
resPst<-NA
for (gg in 1:N){
I<-diag(rep(1,n[gg]))

V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
V.ii[[gg]]<-V[[gg]]%^%(0.5)

resPMp<-matrix(resP[id==gg],ncol=1,nrow=n[gg],byrow=F)
resPMp2<-V.ii.inv[[gg]]%*%resPMp

respermute<-c(respermute,resPMp2)

if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])

resPMpC<-matrix(res.i.c[id==gg],ncol=1,nrow=n[gg],byrow=F)
resPMpC2<-S.i[[gg]]%*%resPMpC
resPMpC2<-resPMpC2

resIst<-c(resIst,resPMpC2)


}

respermute<-respermute[-1]
resIst<-resIst[-1]


resoI2<-resIst[orI]
 t01<- estI

for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
ig<-which(round(t01[orI],10)==round(ii,10))
resoI2[ig]<-sum(resoI2[ig])/length(ig)
}

WI2<-1/sqrt(N )*cumsum(resoI2)

 WsI2 <-list()
 estIm<-list()




iiii=0

while (iiii < M){


if (verbose) print(paste("Iteration: ",iiii+1,sep=""))

if (type=="sign.flip"){

smp<-sample(c(-1,1),size=sum(n),replace=TRUE)

ys<-NA
for (gg in 1:N){
ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
}
ys<-ys[-1] } else {

ys<-NA
for (gg in 1:N){

if (n[gg]==1) smp<-1 else smp<-sample(1:n[gg])
ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute[id==gg])[smp]   )  )
}
ys<-ys[-1]

}

datas<-fit$data
datas[,as.character(fit$call$fixed)[2]]<-ys



fits<-suppressWarnings(update(fit,data=datas))



#sim.proc<-try(get.sim.proc.O.test(fits, residuals=residuals,std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,order.by.original=order.by.original,n=n,N=N,x=x,ZZ=ZZ,id=id, est.original=estI,or.original=orI),silent = TRUE)

#if (class(sim.proc)!="try-error"){
 # iiii=iiii+1
# WsI2[[iiii]]<-sim.proc[[1]]
#estIm[[iiii]]<-sim.proc[[2]]


#}

sim.proc<-get.sim.proc.O.test(fits, residuals=residuals,std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,order.by.original=order.by.original,n=n,N=N,x=x,ZZ=ZZ,id=id, est.original=estI,or.original=orI)
iiii=iiii+1
WsI2[[iiii]]<-sim.proc[[1]]
estIm[[iiii]]<-sim.proc[[2]]


} #end while







res<-list(O=WI2,F=NULL,Om=WsI2,Fm=NULL,Fs=NULL,Fsm=NULL,predO=estI,predOm=estIm,predF=NULL,predFm=NULL,predFs=NULL,predFsm=NULL)
class(res)<-"gofLMM"
res


} #end of function



#' Summary function
#'
#' summary function
#' @param object an object returned by a call to a function
#' @param ... additional arguments passed to from or to other methods


my.summary.gofLMM.testO<-function(object,...){

  O.s<-test.stat.p.val(object$O,object$Om)


  res<-O.s

  rownames(res)<-c(paste("O",rownames(res)[1:2],sep=":"))
  res

}






#' Internal function
#' @keywords internal

get.sim.proc.O.test.2<-function(fit, residuals ,std.type ,use.correction.for.imbalance ,order.by.original,  n,N,x,ZZ,id, orest,ororder ){



  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )

  estI<-fitted(fit,level=1)
  estP<-fitted(fit,level=0)

  orI<-order(estI)



  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()
  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]

  }

  H.i<-solve(H)


  Zb<-matrix(estI-estP,ncol=1)

  res.i.c<-resI

  for (gg in 1:N){
    A<-sigma.est*V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

    B<-Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

    res.i.c[id==gg]<-resI[id==gg]-A%*%ginv(B)%*%Zb[id==gg]

  }


  V.ii.inv<-list()
  V.ii<-list()
  S.i<-list()




  resIst<-NA


  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
    V.ii[[gg]]<-V[[gg]]%^%(0.5)



    if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])

    resPMpC<-matrix(res.i.c[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-S.i[[gg]]%*%resPMpC
    resPMpC2<-resPMpC2

    resIst<-c(resIst,resPMpC2)


  }



  resIst<-resIst[-1]
 if (order.by.original==TRUE) {estI=orest; orI=ororder}

  resoI2<-resIst[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2<-1/sqrt(N )*cumsum(resoI2)

  list(WI2,estI)

}



#' Goodness-of fit test for LMM
#'
#' Goodness-of fit test based on cumulative sum stochastic process for O using a limit expressions for diagonal blocked matrices A and B.

#'
#' @param fit The result of a call to \code{"nlme"}. The model must be fitted with \code{control=lmeControl( returnObject = TRUE)} and \code{keep.data=TRUE}. An error message is returned otherwise. ID variable must be numeric and ordered from 1:N ! Cannot use transformations of the outcome variable directly in the formula i.e. lme(sqrt(y)~x) will return p=1!
#' @param residuals Residuals to be used when constructing the process. Currently implemented only for \code{"individual"} for \emph{individual} residuals.
#' @param std.type Type of standardization to be used for the residuals when constructing the process.
#' Currently implemeneted options are \code{1} and \code{2} for \eqn{S_i=\hat\sigma^{-1/2}I_{n_i}} and \eqn{S_i=\hat{V}_i^{-1/2}}.
#' @param use.correction.for.imbalance Logical. use $n_i^{-1/2} S_i$ when standardizing the residuals. Defaults to \code{FALSE}.
#' @param type How to obtain the processes \eqn{W^m}. Possible values are  \code{"sign.flip"} for the sign-flipping approach and \code{"permutation"} for the permutation approach.
#' @param M Number of random simulations/sign-flipps/permutations. Defaults to \code{100}.
#' @param order.by.original Order the residuals by the original fitted values. Deafults to FALSE.
#' @param verbose Logical. Print the current status of the test. Can slow down the algorithm, but it can make it feel faster. Defaults to \code{FALSE}.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @seealso \code{\link{gof.lmm}}
#' @export

gof.lmm.O.test.2<-function(fit,residuals="individual",std.type=c(1,2),use.correction.for.imbalance=FALSE,type=c("sign.flip","permutation"),M=100,order.by.original=FALSE,verbose=FALSE){

  ####checks, warnings

  if (is.null(fit$data)) stop("Model was fitted with keep.data=FALSE. Use keep.data=TRUE.")

  if (verbose) cat("Using  \"verbose=FALSE \" slows down the algorithm, but it might feel faster. \n")



  ####preliminaries




  id<-fit$data[,names(formula(fit$modelStruct$reStr))]

  N<-length(unique(id))
  n<-table(id)


  id.c<-NA
  for (ii in 1:N){
    id.c<-c(id.c,rep(ii,n[ii]))
  }
  id.c<-id.c[-1]

  if (sum(as.numeric(id)-id.c)!=0) stop("The ID variables needs to be numeric and ordered from 1:N.")

  x<-model.matrix(fit, data=fit$data   )

  ZZ<- model.matrix(formula(fit$modelStruct$reStr)[[1]],data=fit$data)



  ###start gof



  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )

  estI<-fitted(fit,level=1)
  estP<-fitted(fit,level=0)

  orI<-order(estI)
  orP<-order(estP)



  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()
  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]

  }

  H.i<-solve(H)


  Zb<-matrix(estI-estP,ncol=1)

  res.i.c<-resI

    for (gg in 1:N){
          A<-sigma.est*V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

        B<-Z[[gg]]%*%D%*%t(Z[[gg]])%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

        res.i.c[id==gg]<-resI[id==gg]-A%*%ginv(B)%*%Zb[id==gg]

    }






  V.ii.inv<-list()
  V.ii<-list()
  S.i<-list()


  respermute<-NA
  resIst<-NA
  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
    V.ii[[gg]]<-V[[gg]]%^%(0.5)

    resPMp<-matrix(resP[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMp2<-V.ii.inv[[gg]]%*%resPMp

    respermute<-c(respermute,resPMp2)

    if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])

    resPMpC<-matrix(res.i.c[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-S.i[[gg]]%*%resPMpC
    resPMpC2<-resPMpC2

    resIst<-c(resIst,resPMpC2)


  }

  respermute<-respermute[-1]
  resIst<-resIst[-1]


  resoI2<-resIst[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2<-1/sqrt(N )*cumsum(resoI2)

  WsI2 <-list()
  estIm<-list()






  for (iiii in 1:M){

    if (verbose) print(paste("Iteration: ",iiii+1,sep=""))

    if (type=="sign.flip"){

      smp<-sample(c(-1,1),size=sum(n),replace=TRUE)

      ys<-NA
      for (gg in 1:N){
        ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
      }
      ys<-ys[-1] } else {

        ys<-NA
        for (gg in 1:N){

          if (n[gg]==1) smp<-1 else smp<-sample(1:n[gg])
          ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute[id==gg])[smp]   )  )
        }
        ys<-ys[-1]

      }

    datas<-fit$data
    datas[,as.character(fit$call$fixed)[2]]<-ys



    fits<-suppressWarnings(update(fit,data=datas))



    sim.proc<-get.sim.proc.O.test.2(fits, residuals=residuals,std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,order.by.original=order.by.original,n=n,N=N,x=x,ZZ=ZZ,id=id,orest=estI,ororder=orI)

    WsI2[[iiii]]<-sim.proc[[1]]
    estIm[[iiii]]<-sim.proc[[2]]



  } #end for







  res<-list(O=WI2,F=NULL,Om=WsI2,Fm=NULL,Fs=NULL,Fsm=NULL,predO=estI,predOm=estIm,predF=NULL,predFm=NULL,predFs=NULL,predFsm=NULL)
  class(res)<-"gofLMM"
  res


} #end of function





###when you use SF and order by original you do not have to reestimate A and B, this is what I tried here



#' Internal function
#' @keywords internal

get.sim.proc.O.test.type2<-function(fit, residuals ,std.type ,use.correction.for.imbalance ,n,N,x,ZZ,id, est.original,or.original ,A,B){



  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )

  estI<-fitted(fit,level=1)
  estP<-fitted(fit,level=0)

  orI<-order(estI)



  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()
  Xi<-list()
  Zb<-list()
  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
    if (n[gg]!=1) Xi[[gg]]<-x[id==gg,] else Xi[[gg]]<-matrix(x[id==gg,],nrow=1)

  }

  H.i<-solve(H)


  AA<-A
  BB<-B
  Zb<-matrix(estI-estP,ncol=1)


  if (residuals=="individual") res.i.c<-resI-AA%*%ginv(BB)%*%Zb else res.i.c<-resP-(AA+BB)%*%ginv(BB)%*%Zb
  #if (residuals=="individual") res.i.c<-resI-AA%*%ginv(BB)%*%estI else res.i.c<-resP-(AA+BB)%*%ginv(BB)%*%estI




  V.ii.inv<-list()
  V.ii<-list()
  S.i<-list()




  resIst<-NA


  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
    V.ii[[gg]]<-V[[gg]]%^%(0.5)



    if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])

    resPMpC<-matrix(res.i.c[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-S.i[[gg]]%*%resPMpC
    resPMpC2<-resPMpC2

    resIst<-c(resIst,resPMpC2)


  }



  resIst<-resIst[-1]


  estI=est.original
  orI= or.original

  resoI2<-resIst[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2<-1/sqrt(N )*cumsum(resoI2)

  list(WI2,estI)

}




#' Goodness-of fit test for LMM
#'
#' Goodness-of fit test based on cumulative sum stochastic process for O using non-diagonal blocked matrices A and B. I am not reestimating A and B and always ordering by the original fitted values!

#'
#' @param fit The result of a call to \code{"nlme"}. The model must be fitted with \code{control=lmeControl( returnObject = TRUE)} and \code{keep.data=TRUE}. An error message is returned otherwise. ID variable must be numeric and ordered from 1:N ! Cannot use transformations of the outcome variable directly in the formula i.e. lme(sqrt(y)~x) will return p=1!
#' @param residuals Residuals to be used when constructing the process.
#' @param std.type Type of standardization to be used for the residuals when constructing the process.
#' Currently implemeneted options are \code{1} and \code{2} for \eqn{S_i=\hat\sigma^{-1/2}I_{n_i}} and \eqn{S_i=\hat{V}_i^{-1/2}}.
#' @param use.correction.for.imbalance Logical. use \eqn{n_i^{-1/2} S_i} when standardizing the residuals. Defaults to \code{FALSE}.
#' @param type How to obtain the processes \eqn{W^m}. Possible values are  \code{"sign.flip"} for the sign-flipping approach and \code{"permutation"} for the permutation approach.
#' @param M Number of random simulations/sign-flipps/permutations. Defaults to \code{100}.
#' @param verbose Logical. Print the current status of the test. Can slow down the algorithm, but it can make it feel faster. Defaults to \code{FALSE}.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @seealso \code{\link{gof.lmm}}
#' @export


gof.lmm.O.test.type2<-function(fit,residuals="individual",std.type=c(1,2),use.correction.for.imbalance=FALSE,type=c("sign.flip","permutation"),M=100,verbose=FALSE){

  ####checks, warnings

  if (is.null(fit$data)) stop("Model was fitted with keep.data=FALSE. Use keep.data=TRUE.")

  if (verbose) cat("Using  \"verbose=TRUE \" slows down the algorithm, but it might feel faster. \n")



  ####preliminaries




  id<-fit$data[,names(formula(fit$modelStruct$reStr))]

  N<-length(unique(id))
  n<-table(id)


  id.c<-NA
  for (ii in 1:N){
    id.c<-c(id.c,rep(ii,n[ii]))
  }
  id.c<-id.c[-1]

  if (sum(as.numeric(id)-id.c)!=0) stop("The ID variables needs to be numeric and ordered from 1:N.")

  x<-model.matrix(fit, data=fit$data   )

  ZZ<- model.matrix(formula(fit$modelStruct$reStr)[[1]],data=fit$data)



  ###start gof



  resI<-residuals(fit, level = 1  )


  resP<-residuals(fit, level = 0  )

  estI<-fitted(fit,level=1)
  estP<-fitted(fit,level=0)

  orI<-order(estI)
  orP<-order(estP)



  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  D<-getVarCov(fit)

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()
  Xi<-list()
  Zb<-list()
  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
    if (n[gg]!=1) Xi[[gg]]<-x[id==gg,] else Xi[[gg]]<-matrix(x[id==gg,],nrow=1)

  }

  H.i<-solve(H)


  A<-list()
  B<-list()

  res.i.c<-resI

  mm=0
  for (gg in 1:N){
    for (jj in 1:N){
      mm=mm+1


      if (jj==gg){

        zdz<-  Z[[gg]]%*%D%*%t(Z[[gg]])
        cpd<-Xi[[gg]]%*%H.i%*%t(Xi[[gg]])

        #A[[mm]]<-sigma.est*V.i[[gg]]%*%( V[[gg]]- Xi[[gg]]%*%H.i%*%t(Xi[[gg]]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])
        #B[[mm]]<-Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%( V[[gg]]- Xi[[gg]]%*%H.i%*%t(Xi[[gg]]) )%*%V.i[[gg]]%*%Z[[gg]]%*%D%*%t(Z[[gg]])

        vzd<-V.i[[gg]]%*%( V[[gg]]- cpd )%*%V.i[[gg]]%*%zdz

        A[[mm]]<-sigma.est*vzd
        B[[mm]]<-zdz %*%vzd
      } else {

        zdzj<-Z[[jj]]%*%D%*%t(Z[[jj]])
        zdzg<-Z[[gg]]%*%D%*%t(Z[[gg]])
        cpdj<-Xi[[gg]]%*%H.i%*%t(Xi[[jj]])

        #A[[mm]]<--sigma.est*V.i[[gg]]%*%(  Xi[[gg]]%*%H.i%*%t(Xi[[jj]]) )%*%V.i[[jj]]%*%Z[[jj]]%*%D%*%t(Z[[jj]])
        #B[[mm]]<--Z[[gg]]%*%D%*%t(Z[[gg]]) %*%V.i[[gg]]%*%(  Xi[[gg]]%*%H.i%*%t(Xi[[jj]]) )%*%V.i[[jj]]%*%Z[[jj]]%*%D%*%t(Z[[jj]])

        vzdj<-V.i[[gg]]%*%(  cpdj )%*%V.i[[jj]]%*%zdzj

        A[[mm]]<--sigma.est*vzdj
        B[[mm]]<--zdzg %*%vzdj


      }

    }}

  aa<-bb<-list()
  mm=0

  for (gg in 1:N){

    for (jj in 1:N){
      mm=mm+1
      if (jj==1) aa[[gg]]<-A[[mm]] else aa[[gg]]<-cbind(aa[[gg]],A[[mm]])
      if (jj==1) bb[[gg]]<-B[[mm]] else bb[[gg]]<-cbind(bb[[gg]],B[[mm]])
    }}

  for (gg in 1:N){
    if (gg==1) AA<-aa[[gg]] else AA<-rbind(AA,aa[[gg]])
    if (gg==1) BB<-bb[[gg]] else BB<-rbind(BB,bb[[gg]])

  }


  Zb<-matrix(estI-estP,ncol=1)

  if (residuals=="individual") res.i.c<-resI-AA%*%ginv(BB)%*%Zb else res.i.c<-resP-(AA+BB)%*%ginv(BB)%*%Zb
  #if (residuals=="individual") res.i.c<-resI-AA%*%ginv(BB)%*%estI else res.i.c<-resP-(AA+BB)%*%ginv(BB)%*%estI



  V.ii.inv<-list()
  V.ii<-list()
  S.i<-list()


  respermute<-NA
  resIst<-NA
  resPst<-NA
  for (gg in 1:N){
    I<-diag(rep(1,n[gg]))

    V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
    V.ii[[gg]]<-V[[gg]]%^%(0.5)

    resPMp<-matrix(resP[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMp2<-V.ii.inv[[gg]]%*%resPMp

    respermute<-c(respermute,resPMp2)

    if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])

    resPMpC<-matrix(res.i.c[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-S.i[[gg]]%*%resPMpC
    resPMpC2<-resPMpC2

    resIst<-c(resIst,resPMpC2)


  }

  respermute<-respermute[-1]
  resIst<-resIst[-1]


  resoI2<-resIst[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2<-1/sqrt(N )*cumsum(resoI2)

  WsI2 <-list()
  estIm<-list()




  iiii=0

  while (iiii < M){


    if (verbose) print(paste("Iteration: ",iiii,sep=""))

    if (type=="sign.flip"){

      smp<-sample(c(-1,1),size=sum(n),replace=TRUE)

      ys<-NA
      for (gg in 1:N){
        ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute*smp)[id==gg]))
      }
      ys<-ys[-1] } else {

        ys<-NA
        for (gg in 1:N){

          if (n[gg]==1) smp<-1 else smp<-sample(1:n[gg])
          ys<-c(ys,estP[id==gg]+  V.ii[[gg]]%*%( (respermute[id==gg])[smp]   )  )
        }
        ys<-ys[-1]

      }

    datas<-fit$data
    datas[,as.character(fit$call$fixed)[2]]<-ys



    fits<-suppressWarnings(update(fit,data=datas))



    sim.proc<-try(get.sim.proc.O.test.type2(fits, residuals=residuals,std.type=std.type,use.correction.for.imbalance=use.correction.for.imbalance,n=n,N=N,x=x,ZZ=ZZ,id=id, est.original=estI,or.original=orI,A=AA,B=BB),silent = TRUE)

    if (class(sim.proc)!="try-error"){
      iiii=iiii+1
      WsI2[[iiii]]<-sim.proc[[1]]
      estIm[[iiii]]<-sim.proc[[2]]


    }



  } #end while







  res<-list(O=WI2,F=NULL,Om=WsI2,Fm=NULL,Fs=NULL,Fsm=NULL,predO=estI,predOm=estIm,predF=NULL,predFm=NULL,predFs=NULL,predFsm=NULL)
  class(res)<-"gofLMM"
  res


} #end of function






#' Goodness-of fit test for LMM
#'
#' Goodness-of fit test based on cumulative sum stochastic process for O using non-diagonal blocked matrices A and B, simulation approach where refitting is not necessary.

#'
#' @param fit The result of a call to \code{"nlme"}. The model must be fitted with \code{control=lmeControl( returnObject = TRUE)} and \code{keep.data=TRUE}. An error message is returned otherwise. ID variable must be numeric and ordered from 1:N ! Cannot use transformations of the outcome variable directly in the formula i.e. lme(sqrt(y)~x) will return p=1!
#' @param residuals Residuals to be used when constructing the process.
#' @param std.type Type of standardization to be used for the residuals when constructing the process.
#' Currently implemeneted options are \code{1} and \code{2} for \eqn{S_i=\hat\sigma^{-1/2}I_{n_i}} and \eqn{S_i=\hat{V}_i^{-1/2}}.
#' @param use.correction.for.imbalance Logical. use \eqn{n_i^{-1/2} S_i} when standardizing the residuals. Defaults to \code{FALSE}.
#' @param type How to obtain the processes \eqn{W^m}. Possible values are  \code{"sign.flip"} for the sign-flipping random matrix, \code{"normal"} for the standard normal (same for all within cluster), \code{"normal.m"} for the standard normal (different for all within cluster).
#' @param M Number of random simulations/sign-flipps. Defaults to \code{100}.
#' @param verbose Logical. Print the current status of the test. Can slow down the algorithm, but it can make it feel faster. Defaults to \code{FALSE}.
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @seealso \code{\link{gof.lmm}}
#' @export


gof.lmm.O.test.norefit<-function(fit,residuals="individual",std.type=c(1,2),use.correction.for.imbalance=FALSE,type=c("sign.flip","normal","normal.m"),M=100,verbose=FALSE){

  ####checks, warnings

  if (is.null(fit$data)) stop("Model was fitted with keep.data=FALSE. Use keep.data=TRUE.")

  if (verbose) cat("Using  \"verbose=TRUE \" slows down the algorithm, but it might feel faster. \n")



  ####preliminaries

  id<-fit$data[,names(formula(fit$modelStruct$reStr))]

  N<-length(unique(id))
  n<-table(id)


  id.c<-NA
  for (ii in 1:N){
    id.c<-c(id.c,rep(ii,n[ii]))
  }
  id.c<-id.c[-1]

  if (sum(as.numeric(id)-id.c)!=0) stop("The ID variables needs to be numeric and ordered from 1:N.")

  x<-model.matrix(fit, data=fit$data   )
  b<-matrix(c(t(as.matrix(ranef(fit)))),ncol=1)
  resP<-residuals(fit, level = 0  )
  resI<-residuals(fit, level = 1  )


  ZZ<- model.matrix(formula(fit$modelStruct$reStr)[[1]],data=fit$data)

  D<-getVarCov(fit)
  vc<-VarCorr(fit)
  sigma.est<-as.numeric(vc[nrow(vc),1])

  beta.f<-fixef(fit)

  V<-list()
  V.i<-list()
  Z<-list()

  H<-matrix(0,ncol=ncol(x),nrow=ncol(x))
  for (gg in 1:N){
    if (ncol(ZZ)==1) Z[[gg]]<-matrix(ZZ[id==gg,],ncol=1) else Z[[gg]]<-ZZ[id==gg,]
    if (n[gg]==1) Z[[gg]]<-matrix(Z[[gg]],nrow=1)
    I<-diag(rep(1),n[[gg]])
    V[[gg]]<-Z[[gg]]%*%D%*%t(Z[[gg]])+sigma.est*I
    V.i[[gg]]<-V[[gg]]%^%(-1)
    if (n[gg]!=1) H<-H+t(x[id==gg,])%*%V.i[[gg]]%*%x[id==gg,] else H<-H+matrix(x[id==gg,],ncol=1)%*%V.i[[gg]]%*%x[id==gg,]
  }

  H.i<-solve(H)


  ncum<-c(0,cumsum(n))
  Vm<- matrix(0,sum(n),sum(n))
  for (gg in 1:N){
    idr<-idc<- seq(from=1+ncum[gg],to=ncum[gg+1],by=1)
    Vm[idr,idc]<-V[[gg]]
  }

  Vmi<-solve(Vm)

  k<-ncol(D)
  Zm<-matrix(0,sum(n),N*k)

  for (gg in 1:N){
    idr <- seq(from=1+ncum[gg],to=ncum[gg+1],by=1)
    idc<-seq(from= 1+k*(gg-1),to= gg*k  , by=1 )
    Zm[idr,idc]<-Z[[gg]]
  }



  IN<-matrix(0,sum(n),sum(n))
  diag(IN)<-1


  p1<-(Vm-x%*%H.i%*%t(x))%*%(IN- sigma.est*Vmi )
  Am<-sigma.est*Vmi%*%p1
  Bm<-(IN- sigma.est*Vmi )%*%p1

  In<-matrix(0,ncol=N,nrow=N)
  diag(In)<-1
  if (residuals=="individual") J<-sigma.est*Vmi-Am%*%ginv(Bm)%*%Zm%*%(kronecker(In,D))%*%t(Zm)%*%Vmi else J<-IN-(Am+Bm)%*%ginv(Bm)%*%Zm%*%(kronecker(In,D))%*%t(Zm)%*%Vmi


  resT<-J%*%matrix(resP,ncol=1)
  res.i.c<-resT



  V.ii.inv<-list()
  V.ii<-list()
  S.i<-list()


  respermute<-NA
  resIst<-NA
  resPst<-NA
  for (gg in 1:N){
    V.ii.inv[[gg]]<-V[[gg]]%^%(-0.5)
    V.ii[[gg]]<-V[[gg]]%^%(0.5)

    if (std.type==2) S.i[[gg]]<-V.ii.inv[[gg]] else S.i[[gg]]<-  1/sqrt( sigma.est )*diag(rep(1,n[gg]))
    if (use.correction.for.imbalance==TRUE) S.i[[gg]]<-S.i[[gg]]/sqrt(n[gg])

    resPMpC<-matrix(res.i.c[id==gg],ncol=1,nrow=n[gg],byrow=F)
    resPMpC2<-S.i[[gg]]%*%resPMpC
    resPMpC2<-resPMpC2

    resIst<-c(resIst,resPMpC2)


  }

   resIst<-resIst[-1]

  estI<-fitted(fit,level=1)
  estP<-fitted(fit,level=0)

  orI<-order(estI)
  orP<-order(estP)

  resoI2<-resIst[orI]
  t01<- estI

  for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
    ig<-which(round(t01[orI],10)==round(ii,10))
    resoI2[ig]<-sum(resoI2[ig])/length(ig)
  }

  WI2<-1/sqrt(N )*cumsum(resoI2)

  WI2mm<-list()

  ncum<-c(0,cumsum(n))
  Sm<- matrix(0,sum(n),sum(n))
  for (gg in 1:N){
    idr<-idc<- seq(from=1+ncum[gg],to=ncum[gg+1],by=1)
    Sm[idr,idc]<-S.i[[gg]]
  }


  for (kkkk in 1:M){

    if (verbose) print(paste("Iteration: ",kkkk,sep=""))

    respermute<-NA

    for (gg in 1:N){
      I<-diag(rep(1,n[gg]))
      PI<-I
      if (type=="normal") diag(PI)<-rnorm(1)
      if (type=="normal.m") diag(PI)<-rnorm(n[gg])
      if (type=="sign.flip") diag(PI)<-sample(c(-1,1),n[gg],replace=TRUE)

      resPMp<-matrix(resP[id==gg],ncol=1,nrow=n[gg],byrow=F)
      if (type=="normal.m"|type=="sign.flip") resPMp2<-V.ii[[gg]]%*%PI%*%V.ii.inv[[gg]]%*%resPMp else resPMp2<-PI%*%resPMp

      respermute<-c(respermute,resPMp2)



    }

    respermute<-respermute[-1]


    IIp<-matrix(0,ncol(x),1)
    for (gg in 1:N){
      IIp<-IIp+t(x[id==gg,])%*%V.i[[gg]]%*%matrix(respermute[id==gg],ncol=1)

    }




    EMm<-rep(NA,sum(n))
    for (gg in 1:N){
      EMm[id==gg]<-respermute[id==gg]-x[id==gg,]%*%H.i%*%IIp
    }


    resprocsimm<-Sm%*%J%*%EMm
    resoI2mm<-resprocsimm[orI]
    t01<- estI

    for (ii in as.numeric(names(table(t01[orI]))[which(table(t01[orI])>1)])){
      ig<-which(round(t01[orI],10)==round(ii,10))
      resoI2mm[ig]<-sum(resoI2mm[ig])/length(ig)
    }
    WI2mm[[kkkk]]<-1/sqrt(N )*cumsum(resoI2mm)

  }




  res<-list(O=WI2,F=NULL,Om=WI2mm,Fm=NULL,Fs=NULL,Fsm=NULL,predO=estI,predOm=estI,predF=NULL,predFm=NULL,predFs=NULL,predFsm=NULL)
  class(res)<-"gofLMM"
  res


}

