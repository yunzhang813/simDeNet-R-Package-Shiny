\name{sim.expr}
\alias{sim.expr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate expression data
}
\description{
This function simulates expression data for pure and mixed samples, with two cell types (namely, T and N).
}
\usage{
sim.expr(n.samp, mu.T, mu.N, Sigma.T, Sigma.N, prop.T="seq")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.samp}{number of samples to simulate for each cell type.}
  \item{mu.T}{mean vector for cell type T.}
  \item{mu.N}{mean vector for cell type N.}
  \item{Sigma.T}{covariance matrix for cell type T.}
  \item{Sigma.N}{covariance matrix for cell type N.}
  \item{prop.T}{a vector of mixing proportions for cell type T. Options: \code{"seq"} (default), \code{"unif"}, or a customized vector of length \code{n.samp}. If \code{"seq"}, \code{prop.T = seq(0, 1, length=n.samp)}; if \code{"unif"}, \code{prop.T = runif(n.samp, 0, 1)}. And proportions for cell type N are \code{1-prop.T}.} 
}
\details{
The function takes the mean expression (at the log2-transformed scale) and designed covariance structure for two cell types, namely T and N. Mixing proportion of cell type T is required; cell type N takes the complementary proportion. Pure cell-type samples are simulated from multivariate normal (MVN) distribution at the log2-transformed scale. Mixing of cell types is made at the raw scale, i.e. 2^x. The mixed samples are then log2-transformed back.

All MVN simulated expression data are forced to be non-negative by replacing negative values with zero. All simulated samples are quantile normalized by \code{normalize.quantiles()} from the \code{preprocessCore} package.
}
\value{
A list with the following components:
  \item{expr.pure.T}{expression matrix for pure cell type T.}
  \item{expr.pure.N}{expression matrix for pure cell type N.}
  \item{expr.mixed}{expression matrix for mixed cell types.}
}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[mvtnorm]{rmvnorm}} from the \code{mvtnorm} package is used for simulate random MVN samples.

\code{\link[preprocessCore]{normalize.quantiles}} from the \code{preprocessCore} package is used for quantile normalization.
}
\examples{
## load cell type data
data("celltype")
mu.T <- expr[,ctab$Fastq_file_name[which(ctab$X3_letter_code=="ASM")]]
mu.N <- expr[,ctab$Fastq_file_name[which(ctab$X3_letter_code=="AEC")]]

## number of genes
m.gene <- length(mu.T)

## parameters
n.samp <- 5
rho <- c(0.9,0.8,0.7)
block.size <- c(5,10,15)
str.type <- c("interchangeable","decaying","star")

## mixing proportions for cell type T
prop.T <- seq(0, 1, length=n.samp)

## generate covariance matrices
out.Sigma.T <- generate.Sigma(m.gene=m.gene, dd=rep(1,m.gene), rho=rho, block.size=block.size, str.type=str.type)
Sigma.T <- out.Sigma.T$Sigma
Sigma.N <- diag(m.gene) # identity matrix

## simulate pure and mixed expression samples
out.expr <- sim.expr(n.samp, mu.T, mu.N, Sigma.T, Sigma.N, prop.T)

## If mixing proportions for cell type T are all 0's
prop0 <- rep(0, n.samp)
## mixed samples should be the same as the pure samples for cell type N
out.expr0 <- sim.expr(n.samp, mu.T, mu.N, Sigma.T, Sigma.N, prop0)
identical(as.numeric(out.expr0$expr.mixed), as.numeric(out.expr0$expr.pure.N))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }% use one of  RShowDoc("KEYWORDS")

