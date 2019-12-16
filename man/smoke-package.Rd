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
\keyword{ package }
\examples{
\donttest{data.raw=readTraces("traces.csv")}
\donttest{data.scaled=scaLoad(data.raw)}
\donttest{data.doubleRefed=doubleRef(data.scaled)}
\donttest{rate0=estKinetics(data.doubleRefed)}
\donttest{kinetics=fitKinetics(data.doubleRefed)}
\donttest{kinetics=fitKinetics(data.doubleRefed,plotResidual=TRUE)}
}
