\name{doubleBlank}
\alias{doubleBlank}
\title{
Substract the double-references 
}
\description{
The BLI experiment for small-molecules is designed to have double references: reference biosensor and reference sample, i.e. b
}
\usage{
doubleBlank(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
}
\details{
The signal in small-molecule BLI experiment is small in relation to noise. So the experiment design includes double references: reference biosensor and reference sample (no ligand). With subtraction of both, we expect to remove both the background signal and the non-specific binding. \cr
The \code{doubleBlank} assumes that \code{traces} in the input "Bli" object contains two consecutive sets of columns with the second set as the refence biosensor, and the reference samples as the last column in each set. This means that the traces are "paired" and traces number is even. If not, an error message will be raised: "the double-referenced traces should be paired".\cr
After \code{doubleBlank}, the number of traces should be equal to the number of ligand concentration. If not, a warning will be raised: "ligand concentrations mis-match trace number".\cr
Also, after successful \code{doubleBlank}, the \code{status} of "doubleBlank" will be set to "TRUE". This prevent the user to run \code{doubleBlank} twice in accident.
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
\donttest{load("bli.RData")}
\donttest{bli1 <- alignLoad(bli, 180, 780)}
\donttest{bli2 <- doubleBlank(bli1)}
}
\keyword{ manip }
