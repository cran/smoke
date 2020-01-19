\name{estimate}
\alias{estimate}
\title{
Estimate initial kinetic parameters
}
\description{
Estimate the initial kinetic parameters, \code{kOn0} and \code{kOff0}, and save in the returned "Bli" object. These parameters can serve as the starting value in the next step, \code{fitKinetics}, to calculate the kinetics binding model of the small-molecule under investigation.
}
\usage{
estimate(obj)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
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
\note{
\code{kOn0} and \code{kOff0} can also be input manually with the setter functions, \code{\link{kOn0<-}} and \code{\link{kOff0<-}}.
}
\examples{
\donttest{load("bli.RData")}
\donttest{bli1 <- alignLoad(bli, 180, 780)}
\donttest{bli2 <- doubleBlank(bli1)}
\donttest{bli3 <- baseline(bli2, 1080, 1260)}
\donttest{bli4 <- estimate(bli3)}
}
\keyword{ manip }
