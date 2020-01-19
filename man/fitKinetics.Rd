\name{fitKinetics}
\alias{fitKinetics}
\title{
Fit binding kinetics
}
\description{
Fit the BLI traces with a single association-then-dissociation equation.
}
\usage{
fitKinetics(obj)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
}
\details{
We fit the whole set of BLI traces with a single association-then-dissociation equation \insertCite{pza}{smoke}. The fitting model will be stored in the \code{kinetics} slot, and the \code{status} of "fitKinetics" will be set to "TRUE".
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
\seealso{
  \code{\link{kinetics}}, \code{\link{plotKinetics}}, \code{\link{plotResiduals}}
}
\examples{
\donttest{load("bli.RData")}
\donttest{bli1 <- alignLoad(bli, 180, 780)}
\donttest{bli2 <- doubleBlank(bli1)}
\donttest{bli3 <- baseline(bli2, 1080, 1260)}
\donttest{bli4 <- estimate(bli3)}
\donttest{bli5 <- fitKinetics(bli4)}
}
\keyword{ manip } 
