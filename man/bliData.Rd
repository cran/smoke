\name{bli}
\alias{bli}
\docType{data}
\title{
BLI data example
}
\description{
This is an unprocessed BLI data saved in a "BLI" object \insertCite{pza}{smoke}.
}
\usage{data("bli")}
\format{
  A S4 "Bli" class with 7 slots.
  \describe{
    \item{\code{status}}{a logical vector indicating this is a raw dataset}
    \item{\code{traces}}{a data.frame with 16 BLI traces at 12290 time points}
    \item{\code{ligand}}{a numeric vector of ligand concentrations}
    \item{\code{tExp}}{a numeric vector of association and dissociation times}
    \item{\code{kOn0}}{a numeric of initial on-rate; empty} 
    \item{\code{kOff0}}{a numeric of initial off-rate; empty}
    \item{\code{kinetics}}{a \code{nls} object; empty} 
  }
}
\source{
\insertRef{pza}{smoke}
}
\examples{
data(bli)
}
\keyword{datasets}

