\name{plotResiduals}
\alias{plotResiduals}
\title{
Plot residuals from kinetics model fitting
}
\description{
Plot the residuals from the kinetics-model fitting
}
\usage{
plotResiduals(obj, ...)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
  \item{\dots}{
Arguments to be passed to methods, such as graphical parameters
}
}
\references{
\insertRef{pza}{smoke}
}
\author{
\packageAuthor{smoke}

Maintainer: \packageMaintainer{smoke}
}
\note{
This method only works after \code{fitKinetics}.
}
\seealso{
  \code{\link{fitKinetics}}, \code{\link{plotKinetics}}
}
\examples{
\donttest{load("bli.RData")}
\donttest{bli1 <- alignLoad(bli, 180, 780)}
\donttest{bli2 <- doubleBlank(bli1)}
\donttest{bli3 <- baseline(bli2, 1080, 1260)}
\donttest{bli4 <- estimate(bli3)}
\donttest{bli5 <- fitKinetics(bli4)}
\donttest{plotResiduals(bli5)}
}
\keyword{ hplot }
