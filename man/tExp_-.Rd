\name{tExp<-}
\alias{tExp<-}
\title{
Setter for BLI times \code{tExp}
}
\description{
The setter function for BLI times \code{tExp}, the start time of association and dissociation steps.
}
\usage{
tExp(obj) <- value
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
  \item{value}{
"numeric" vector for the start time of association and dissociation steps
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
\examples{
bli <- Bli()
tExp(bli) <- c(1260,1860)
}
\keyword{ manip }
