\name{alignLoad}
\alias{alignLoad}
\title{
Align BLI traces with loading step
}
\description{
Align BLI traces with loading step with the start and end times.
}
\usage{
alignLoad(obj, loadStart, loadEnd)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
  \item{loadStart}{
  start time of the loading step
}
  \item{loadEnd}{
  end time of the loading step
}
}
\details{
 This is the optional first step in BLI data processing. The \code{traces} in the input "Bli" object is assumed to be paired original dataset. The output "Bli" object contains the \code{traces} aligned, and the \code{status} of "alignLoad" set to "TRUE".
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
}
\keyword{ manip }
