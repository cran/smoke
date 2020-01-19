\name{kinetics}
\alias{kinetics}
\title{
Output kinetics model
}
\description{
Output kinetics model with matrix containing KD, rMax, kOn, kOff, and the relevant statistics.
}
\usage{
kinetics(obj)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
}
\value{
 "matrix" with kinetic parameters
}
\references{
\insertRef{pza}{smoke}
}
\author{
\packageAuthor{smoke}

Maintainer: \packageMaintainer{smoke}
}
\seealso{
 \code{\link{fitKinetics}}
}
\examples{
\donttest{load("bli.RData")}
\donttest{bli1 <- alignLoad(bli, 180, 780)}
\donttest{bli2 <- doubleBlank(bli1)}
\donttest{bli3 <- baseline(bli2, 1080, 1260)}
\donttest{bli4 <- estimate(bli3)}
\donttest{bli5 <- fitKinetics(bli4)}
\donttest{kinetics(bli5)}
}
\keyword{ manip }
