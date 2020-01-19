\name{ligand}
\alias{ligand}
\title{
Getter for "Bli" slot of ligand concentration
}
\description{
The getter function for "Bli" slot of ligand concentration. If the slot is empty, there will be an error message: "please input ligand concentration".
}
\usage{
ligand(obj)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
}
\value{
a numeric vector for the ligand concentration in BLI experiment
}
\author{
\packageAuthor{smoke}

Maintainer: \packageMaintainer{smoke}
}
\seealso{
 \code{\link{ligand<-}}
}
\examples{
\donttest{load("bli.RData")}
\donttest{conc <- ligand(bli)}
}
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
