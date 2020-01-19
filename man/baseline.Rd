\name{baseline}
\alias{baseline}
\title{
Align BLI traces with baseline step
}
\description{
Align BLI traces with the baseline step right before the association step. The start and end times of baseline are input arguments.
}
\usage{
baseline(obj, tStart, tEnd)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
  \item{tStart}{
  start time of the baseline
}
  \item{tEnd}{
  end time of the baseline
}
}
\details{
 This step follows \code{doubleBlank} step. The output "Bli" object contains the \code{traces} aligned with baselines, and the \code{status} of "baseline" set to "TRUE". Also the traces number will be checked. If it does not equal the number of ligand concentration, a warning will be raised: "ligand concentrations mis-match trace number".
}
\value{
 S4 object of "Bli" class
}
\references{
\insertRef{pza}{smoke}
}
\author{
\packageAuthor{smoke}

Maintainer: \packageMaintainer{smoke}
}
\examples{
\donttest{load("bli.RData")}
\donttest{bli1 <- alignLoad(bli, 180, 780)}
\donttest{bli2 <- doubleBlank(bli1)}
\donttest{bli3 <- baseline(bli2, 1080, 1260)}
}
\keyword{ manip }
