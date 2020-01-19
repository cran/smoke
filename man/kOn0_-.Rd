\name{kOn0<-}
\alias{kOn0<-}
\title{
Setter for "Bli" slot of \code{kOn0}
}
\description{
The setter function for "Bli" slot of \code{kOn0}, initial kOn rate.
}
\usage{
kOn0(obj) <- value
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
  \item{value}{
  initial kOn rate
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
It is preferred to \code{estimate} initial kOn rate.
}
\seealso{
 \code{\link{estimate}}
}
\examples{
bli <- Bli()
kOn0(bli) <- 0.1
}
\keyword{ manip }
