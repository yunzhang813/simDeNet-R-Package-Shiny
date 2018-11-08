\name{deconv}
\alias{deconv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deconvolution
}
\description{
An integrated function for deconvolution of mixed cell-type samples, based on the \code{ISOpureR} package.
}
\usage{
deconv(mixed, ref, seed=123, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mixed}{data matrix (at the log2-transformed scale) for mixed samples, with genes in rows and samples in columns.}
  \item{ref}{data matrix (at the log2-transformed scale) for reference samples, with genes in rows and samples in columns.}
  \item{seed}{random seed for reproducible result.}
  \item{...}{arguments passed to \code{ISOpure.step1.CPE()}.}
}
\value{
A list of the following components:
  \item{expr.deconv}{estimated expression profiles (at the log2-transformed scale) for the targeted pure samples after deconvolution.}
  \item{est.prop}{estimated mixing proportions for the targeted pure samples.}

}
\references{
  G Quon, S Haider, AG Deshwar, A Cui, PC Boutros, QD Morris. \emph{Computational purification of 
  individual tumor gene expression profiles}. Genome Medicine (2013) 5:29, 
  \url{http://genomemedicine.com/content/5/3/29}.
  
  G Quon, QD Morris. \emph{ISOLATE: a computational strategy for identifying the primary origin of
  cancers using high-thoroughput sequencing}. Bioinformatics 2009, 25:2882-2889 
  \url{http://bioinformatics.oxfordjournals.org/content/25/21/2882}.
}
\note{
The original functions in the \code{ISOpureR} package input data at the raw scale (i.e. 2^x) and return the estimated profiles also at the raw scale. Our integrated function takes care the data transformation internally, so that it inputs and returns data at the log2-transformed scale, keeping consistancy with other functions in our package.
}

\section{Warning }{
Use \code{warnings()} to see warnings. Warnings are expected in \code{ISOpureR} from the optimization calculations.
}

\seealso{
\code{\link[ISOpureR]{ISOpure.step1.CPE}} and \code{\link[ISOpureR]{ISOpure.step2.PPE}} from the \code{ISOpureR} package.
}
\examples{
## Firstly, let's generate some data.
set.seed(999)
data("celltype")
mu.T <- expr[,ctab$Fastq_file_name[which(ctab$X3_letter_code=="ASM")]]
mu.N <- expr[,ctab$Fastq_file_name[which(ctab$X3_letter_code=="AEC")]]
## number of samples to simulate
n.samp <- 5
## parameters for correlation design of cell type T
rho <- c(0.9,0.8,0.7)
block.size <- c(5,10,15)
str.type <- c("interchangeable","decaying","star")
## one-step simulation
out.oneStepSim <- oneStepSim(n.samp, mu.T, mu.N, rho=rho, block.size=block.size, str.type=str.type)

## Deconvolution. Pure samples for cell type N is used as reference.
\dontrun{
out.deconv <- deconv(mixed=out.oneStepSim$expr.mixed, ref=out.oneStepSim$expr.pure.N)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }% use one of  RShowDoc("KEYWORDS")
