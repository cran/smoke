\name{kOff0<-}
\alias{kOff0<-}
\title{
Setter for "Bli" slot of \code{kOff0}
}
\description{
The setter function for "Bli" slot of \code{kOff0}, initial kOff rate
}
\usage{
kOff0(obj) <- value
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
  \item{value}{
  initial kOff rate
}
}
\value{
  S4 object of "Bli" class 
}
\author{
\packageAuthor{smoke}

Maintainer: \packageMaintainer{smoke}
}
\note{
It is preferred to \code{estimate} initial kOff rate.
}
\seealso{
 \code{\link{estimate}}
}
\examples{
bli <- Bli()
kOff0(bli) <- 0.1
}
\keyword{ manip }
