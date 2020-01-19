\name{traces}
\alias{traces}
\title{
Getter of BLI \code{traces}
}
\description{
The getter function of BLI \code{traces}. The first column of this data.frame is the time in second; the following columns are the BLI traces. 
}
\usage{
traces(obj)
}
\arguments{
  \item{obj}{
  \code{obj} of "Bli" class
}
}
\value{
"data.frame" of BLI \code{traces}
}
\author{
\packageAuthor{smoke}

Maintainer: \packageMaintainer{smoke}
}
\examples{
\donttest{load("bli.RData")}
\donttest{traces(bli)}
}
\keyword{ manip }

