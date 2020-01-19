\name{status}
\alias{status}
\title{
Getter for data processing \code{status}
}
\description{
The getter function for the data processing \code{status}.
}
\usage{
status(obj)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
}
\value{
"logical" vector as the indicator of the status of data analysis
}
\author{
\packageAuthor{smoke}

Maintainer: \packageMaintainer{smoke}
}
\note{
There is no Setter for \code{status} by design.\cr
Because \code{status} is modified automatically during data processing, there is no need for the user to change the \code{status} manually.
}
\seealso{
  \code{\link{alignLoad}}, \code{\link{doubleBlank}}, \code{\link{baseline}}, \code{\link{estimate}}, \code{\link{fitKinetics}}
}
\examples{
bli <- Bli()
status(bli)
}
\keyword{ manip }
