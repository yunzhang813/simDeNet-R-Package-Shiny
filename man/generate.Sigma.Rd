\name{generate.Sigma}
\alias{generate.Sigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate covariance matrix
}
\description{
This function generates a matrix with block-wise customized structure.
}
\usage{
generate.Sigma(m.gene, block.size, rho, dd=NULL, str.type="interchangeable")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m.gene}{number of genes in the covariance matrix.}
  \item{block.size}{vector of sizes (i.e. number of genes) for each block. Size can not be 0, at least 1.}
    \item{rho}{vector of correlation coefficients (i.e. edge strength) in each block. The vector should have the same length as \code{block.size}. Positive correlations are suggested; use with caution for negative correlations.}
  \item{dd}{diagnal entries of the weight matrix. It should have length equal to \code{m.gene}. If \code{NULL}, all 1's will be assigned.}
  \item{str.type}{types of block structure. Options: \code{"interchangeable"} (default) assumes an interchangeable correlation in the block; \code{"decaying"} is an AR(1) model; \code{"star"} gives a star-shape structure. It can also be a vector of mixed types, but should have the same length as \code{rho}.}
}

\value{
A list with the following components:
  \item{Sigma}{a covariance matrix with customized structure.}
  \item{true.str}{true structure in generated \code{Sigma}. A matrix with 0 (no structure) and 1 (an edge).}
}

\examples{
## a correlation matrix with three designed blocks
Sigma.str1 <- generate.Sigma(m.gene=35, block.size=c(5,10,15), rho=c(0.9,0.8,0.7),
                             str.type=c("interchangeable","decaying","star"))
                             
## no structure, i.e. an identity matrix
Sigma.str2 <- generate.Sigma(m.gene=10, block.size=1, rho=0)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }% use one of  RShowDoc("KEYWORDS")
