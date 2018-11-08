\name{celltype}
\alias{celltype}
\docType{data}
\title{
The cell type data
}
\description{
Sample dataset of pure cell type microRNA abundance
}
\usage{data("celltype")}
\format{
  A data frame with 382 microRNAs and 32 cell-type specific samples.
  \describe{
    \item{\code{expr}}{data matrix with 382 rows (miRNAs) and 32 columns (samples).}
    \item{\code{ctab}}{column information of \code{expr}.}
  }
}
\examples{
data(celltype)
}
\keyword{datasets}
