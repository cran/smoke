\name{smoke-package}
\alias{smoke-package}
\alias{smoke}
\docType{package}
\title{
\packageTitle{smoke}
}
\description{
\packageDescription{smoke}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{smoke}
\packageIndices{smoke}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
\packageAuthor{smoke}

Maintainer: \packageMaintainer{smoke}
}
\references{
\insertRef{pza}{smoke}
}
\keyword{ package }
\examples{
\donttest{load("bli.RData")}
\donttest{bli1 <- alignLoad(bli, 180, 780)}
\donttest{bli2 <- doubleBlank(bli1)}
\donttest{bli3 <- baseline(bli2, 1080, 1260)}
\donttest{bli4 <- estimate(bli3)}
\donttest{bli5 <- fitKinetics(bli4)}
\donttest{plotTraces(bli2)}
\donttest{plotKinetics(bli5)}
\donttest{plotResiduals(bli5)}
}

