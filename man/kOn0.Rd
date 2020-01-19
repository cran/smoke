\name{kOn0}
\alias{kOn0}
\title{
Getter for "Bli" slot of \code{kOn0}
}
\description{
The getter function for "Bli" slot of \code{kOn0}, initial kOn rate.
}
\usage{
kOn0(obj)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
}
\value{
numeric of \code{kOn0}
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
\donttest{kOn0(bli4)}
}
\keyword{ manip }
