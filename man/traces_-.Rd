\name{traces<-}
\alias{traces<-}
\title{
Setter of BLI \code{traces}
}
\description{
The setter function of BLI \code{traces}. The input "data.frame" contains the times in the first column and the BLI traces in the following columns.
}
\usage{
traces(obj) <- value
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
  \item{value}{
"data.frame" of BLI \code{traces}
}
}
\details{
There may be two kinds of input "data.frame":\cr
The "original" data with times in the first column, and two consecutive sets of BLI traces in the remaining columns. The second set are the traces from the reference biosensor. In both sets, the last trace is from the reference sample (no ligand).\cr
The user may input the traces of which the double references have been subtracted. In this case, the first column is still the time. The rest of the columns contain the traces with the exact number of the length of ligand concentration.
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
\seealso{
 \code{\link{doubleBlank}}
}
\examples{
\donttest{bli <- Bli()}
\donttest{# traces(bli) <- read.csv("traces.csv")}
}
\keyword{ manip }
