\name{eval.ROC}
\alias{eval.ROC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate ROC and AUC
}
\description{
This function performs ROC analysis for network estimation based on the \code{ROCR} package.
}
\usage{
eval.ROC(est.str, true.str, plot.ROC=TRUE, show.AUC=TRUE, zoom=FALSE, lightPDF=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{est.str}{symmetric matrix or array of symmetric matrices of estimated network structure, with elements taking values form 0 to 1 which measures the possibility of being an edge. It should have the same dimension as \code{true.str}.}
  \item{true.str}{symmetric matrix or array of symmetric matrices of true structure, usually with binary elements of 1 = edge and 0 = no edge. It should have the same dimension as \code{est.str}.}
  \item{plot.ROC}{logical value indicating if to plot the ROC curve. Default: \code{TRUE}.}
  \item{show.AUC}{logical value indicating if to show AUC on the plot. Default: \code{TRUE}.}
  \item{zoom}{logical value indicating if to zoom-in the plot. Default: \code{FALSE}.}
  \item{lightPDF}{logical value indicating if to subsample points for plotting in a PDF. Default: \code{FALSE}.}
  \item{...}{arguments passed to \code{plot()}.}
}

\value{
A list with the following components:
  \item{perf}{object of \code{performence-class}.}
  \item{AUC}{numeric value of the AUC statistic.}
  \item{xval}{maximum value plotted in the x-axis.}
  \item{yval}{maximum value plotted in the y-axis.}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


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

## If we consider correlation network
est.str.pure <- abs(cor(t(out.oneStepSim$expr.pure.T)))
est.str.mixed <- abs(cor(t(out.oneStepSim$expr.mixed)))
## the true structure in cell type T should be
true.str <- out.oneStepSim$true.str.T

## ROC analysis
out.ROC.pure <- eval.ROC(est.str=est.str.pure, true.str=true.str, plot.ROC=TRUE, main="Pure samples")
out.ROC.pure$AUC
out.ROC.mixed <- eval.ROC(est.str=est.str.mixed, true.str=true.str, plot.ROC=TRUE, main="Mixed samples")
out.ROC.mixed$AUC
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dplot }% use one of  RShowDoc("KEYWORDS")
