\name{tExp}
\alias{tExp}
\title{
Getter for BLI times \code{tExp}
}
\description{
The getter function for BLI times \code{tExp}, the start time of association and dissociation steps.
}
\usage{
tExp(obj)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
}
\value{
"numeric" vector for the start time of association and dissociation steps
}
\references{
\insertRef{pza}{smoke}
}
\author{
\packageAuthor{smoke}

Maintainer: \packageMaintainer{smoke}
}
\examples{
\donttest{load("bli.RData")}
\donttest{tExp(bli)}
}
\keyword{ manip }
