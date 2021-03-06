\name{lhreg}
\alias{lhreg}
\alias{logLik.lhreg}
\alias{summary.lhreg}
\alias{simulate.lhreg}
\alias{profile_lambda1}
\alias{loo1}
\alias{loo2}
\alias{parametric_bootstrap}
\alias{pred_int}
\title{
Phylogeny and Species Trait Effects on Detectability
}
\description{
Functions used in the manuscript.
}
\usage{
lhreg(Y, X, SE, V, init=NULL, lambda=NA, method="Nelder-Mead",
    hessian=FALSE, DElimit=10, eval=FALSE)

\method{logLik}{lhreg}(object, ...)
\method{summary}{lhreg}(object, ...)
\method{simulate}{lhreg}(object, nsim = 1, seed = NULL,
    lambda = NA, obs_error = FALSE, ...)
profile_lambda1(object, value, ...)

loo1(i, object, return_coefs=TRUE)
loo2(i, object, return_coefs=TRUE, method=NULL)

parametric_bootstrap(object, nsim=1, seed = NULL, method, cl=NULL, ...)
pred_int(object, boot, cl=NULL)
}
\arguments{
  \item{Y}{response vector.}
  \item{X}{model matrix for the mean.}
  \item{SE}{standard error estimate (observation error) for the response.}
  \item{V}{correlation matrix.}
  \item{init}{initial values or \code{NULL}.}
  \item{lambda}{phylogeny strength, non-negative,
    \code{NA} means it is not fixed value but estimated.}
  \item{method}{method argument accepted by \code{\link{optim}},
    or \code{"DE"} when DEoptim is called. Method for \code{loo2} is
    same as for the input \code{object$method} by default (\code{method=NULL}).}
  \item{hessian}{logical, if the Hessian needs to be estimated at MLE.}
  \item{DElimit}{limit for DEoptim search (used as [-DElimit, +DElimit])
    when \code{method = "DE"}.}
  \item{eval}{logical, the negative log-likelihood is evaluated
    at \code{init} given the data arguments without optimization.}
  \item{object}{a fitted object returned by \code{lhreg}.}
  \item{value}{fixed value for \code{lambda}
    in profile likelihood calculations.}
  \item{i}{index of observations to drop for cross-validation.}
  \item{nsim}{number of response vectors to simulate. Defaults to 1.}
  \item{seed}{an object specifying if and how the random number generator should
    be initialized as described in \code{\link{simulate}}.}
  \item{obs_error}{logical, if observation error is to be taken into account.}
  \item{cl}{number of parallel processes or cluser object.}
  \item{return_coefs}{logical, if (re)estimated coefficients are to be returned.}
  \item{boot}{an object with \code{parametric_bootstrap} results.}
  \item{\dots}{other arguments passed to underlying functions.}
}
\details{
See Examples and Vignettes for details.
}
\value{
\code{lhreg} returns an object of lhreg, that is a list.

The \code{summary} method returns a summary for the input object.
The \code{logLik} method returns the log-likelihood.
The \code{simulate} method returns the random deviates
under a multivariate normal model.

\code{profile_lambda1} returns log-likelihood based on fixed lambda profile
likelihood.

\code{loo1} returns the observed value for the held-out data point and the
corresponding prediction based on multiple linear regression.
Also returns coefficients based on the training data when
\code{return_coefs=TRUE}.

\code{loo2} returns the observed value for the held-out data point and the
corresponding prediction based on correlated mixed-effects model.
Also returns coefficients based on the training data when
\code{return_coefs=TRUE}.

\code{parametric_bootstrap} uses the \code{simulate} method to simulate
observations from a Multivariate Normal distribution according to the
input object (without the observation error) to refit the model and returns
simulated values and estimates.

\code{pred_int} calculates the prediction interval for an observation
conditional on the other species and the known tree (this one and the other
species included), and returns the bootstrap distribution of the prediction
that can be used to calculate quantile based prediction intervals.
}
\author{
Peter Solymos <solymos@ualberta.ca>
}
\seealso{
\code{\link{optim}}, \code{\link{DEoptim}}
}
\examples{
## see examples in the vignette
\dontrun{
vignette(topic = "lhreg", package = "lhreg")
}
}
\keyword{ models }
