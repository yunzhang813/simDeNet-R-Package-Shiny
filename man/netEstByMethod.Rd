\name{netEstByMethod}
\alias{netEstByMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Network estimation function
}
\description{
This function takes simulated data and estimates network based on WGCNA (correlation-based), ARACNE (mutual information-based) or Bayesian algorithms.
}
\usage{
netEstByMethod(data.list, method, softPower=1)
}
\arguments{
  \item{data.list}{a list of simulated data, e.g. \code{list(data.pure, data.mixed, data.deconvoluted)}. The output of this function will be a list of estimated networks of the same order.}
  \item{method}{method used to estimate network. Options: \code{WGCNA} is correlation-based, \code{ARACNE} is mutual information-base, \code{Bayesian} is Bayesian network, and \code{GLASSO}.}
  \item{softPower}{soft power used in WGCNA. Default: \code{softPower=1}.}
}
\details{
This function uses \code{adjacency()} for \code{method="WGCNA"}, \code{minet()} for \code{method="ARACNE"},  \code{hc(),bn.fit()} for \code{method="Bayesian"}, and simply \code{var()} for \code{method="GLASSO"} for network estimation.
}
\value{
This function returns a list of estimated networks, coresponding to the input data list.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
\dontrun{
## simulated data
data.list <- list(out.oneStepSim$expr.pure.T, out.oneStepSim$expr.mixed, out.deconv$expr.deconv)
names(data.list) <- c("pure", "mixed", "deconvoluted")

## one step analysis
net.WGCNA <- netEstByMethod(data.list, method="WGCNA")
net.ARACNE <- netEstByMethod(data.list, method="ARACNE")
net.Bayesian <- netEstByMethod(data.list, method="Bayesian")
net.GLASSO <- netEstByMethod(data.list, method="GLASSO")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }% use one of  RShowDoc("KEYWORDS")

