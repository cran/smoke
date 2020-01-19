\name{kOff0}
\alias{kOff0}
\title{
Getter for "Bli" slot of \code{kOff0}
}
\description{
The getter function for "Bli" slot of \code{kOff0}, initial kOff rate.
}
\usage{
kOff0(obj)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
}
\value{
numeric of \code{kOff0}
}
\author{
\packageAuthor{smoke}

Maintainer: \packageMaintainer{smoke}
}
\seealso{
 \code{\link{estimate}}
}
\examples{
\donttest{load("bli.RData")}
\donttest{bli1 <- alignLoad(bli, 180, 780)}
\donttest{bli2 <- doubleBlank(bli1)}
\donttest{bli3 <- baseline(bli2, 1080, 1260)}
\donttest{bli4 <- estimate(bli3)}
\donttest{kOff0(bli4)}
}
\keyword{ manip }
