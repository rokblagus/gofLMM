% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_gofLMM.R
\name{sim.data.cluster}
\alias{sim.data.cluster}
\title{Goodness-of fit test for LMM, function to simulate data}
\usage{
sim.data.cluster(
  N,
  n,
  betas,
  norm.eps,
  var.eps = NULL,
  shape = NULL,
  scale = NULL,
  norm.re.intercept,
  var.re.intercept = NULL,
  shape.re.intercept = NULL,
  scale.re.intercept = NULL,
  sim.re.slope,
  norm.re.slope = NULL,
  var.re.slope = NULL,
  shape.re.slope = NULL,
  scale.re.slope = NULL,
  sim.x2.qdr = FALSE,
  b.qdr = NULL
)
}
\arguments{
\item{N}{number of clusters}

\item{n}{number of subjects per cluster (the same for all clusters)}

\item{betas}{Vector of true regression coefficients for the fixed effects}

\item{norm.eps}{Logical, if TRUE the errors are simulated from a (zero mean) normal distribution with variance \code{var.eps}, otherwise from a (centered) gama with parameters \code{shape} and \code{scale}}

\item{var.eps}{see above}

\item{shape}{see above}

\item{scale}{see above}

\item{norm.re.intercept}{Logical, if TRUE the random intercepts are simulated from a (zero mean) normal distribution with variance \code{var.re.intercept}, otherwise from a (centered) gama with parameters \code{shape.re.intercept} and \code{scale.re.intercept}}

\item{var.re.intercept}{see above}

\item{shape.re.intercept}{see above}

\item{scale.re.intercept}{see above}

\item{sim.re.slope}{Logical. If TRUE random slopes are simulated.}

\item{norm.re.slope}{Logical, if TRUE the random slopes are simulated from a (zero mean) normal distribution with variance \code{var.re.slope}, otherwise from a (centered) gama with parameters \code{shape.re.slope} and \code{scale.re.slope}}

\item{var.re.slope}{see above}

\item{shape.re.slope}{see above}

\item{scale.re.slope}{see above}

\item{sim.x2.qdr}{Logical. If TRUE the square of X2 is included in the true (correct) fixed effects design matrix.}

\item{b.qdr}{True beta coefficient associated with the square of X2}
}
\description{
This function can be used to simulate (balanced) cluster data as used in the simulation study of Peterlin et al. See the paper for details.
}
\seealso{
\code{\link{gof.lmm}}
}
\author{
Rok Blagus, \email{rok.blagus@mf.uni-lj.si}
}
