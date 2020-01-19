\name{plotTraces}
\alias{plotTraces}
\title{
Plot BLI traces
}
\description{
Plot the BLI traces in the "Bli" object.
}
\usage{
plotTraces(obj, ...)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
  \item{\dots}{
Arguments to be passed to methods, such as graphical parameters
}
}
\author{
\packageAuthor{smoke}

Maintainer: \packageMaintainer{smoke}
}
\seealso{
  \code{\link{plotKinetics}}, \code{\link{plotResiduals}}
}
\examples{
\donttest{load("bli.RData")}
\donttest{plotTraces(bli)}
}
\keyword{ hplot }
