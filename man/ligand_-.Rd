\name{ligand<-}
\alias{ligand<-}
\title{
Setter for "Bli" slot of ligand concentration
}
\description{
The setter function for "Bli" slot of ligand concentration in BLI experiment.
}
\usage{
ligand(obj) <- value
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
  \item{value}{
  \code{value} numeric vector of ligand concentrations
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
The ligand concentrations do not include the blank control (0).
}


\seealso{
 \code{\link{ligand}}
}
\examples{
\donttest{bli1 <- Bli()}
\donttest{ligand(bli1) <- 16/2^(0:6)}
}
\keyword{ manip }
