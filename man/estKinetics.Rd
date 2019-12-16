\name{estKinetics}
\alias{estKinetics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
estimate initial values of kOn and kOff
}
\description{
A robust and rough way to estimate initial values of kON and kOff from data only. There is no manual input for the initial values
}
\usage{
estKinetics(data, tD0 = 600, lig = c(16/2^(0:6)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the BLI traces with double-references corrected and aligned at start-time of association
}
  \item{tD0}{
start-Time of dissociation step
}
  \item{lig}{
concentrations of small-molecule analyte
}
}
\details{
estimates of the association and dissociation with two steps: from single BLI trace, first log transform the first 50s of the association/dissociation step, linear regression for the slope; secondly use that slope as start value for non-linear fitting for the association/dissociation rate. The mean values from all traces are returned. The fitting curves are ploted. 
}
\value{
numeric vector of estimated values of kOff and kOn
}
\author{
Qingan Sun, Xiaojun Li, James C Sacchettini
}
\note{
Due to the small SNR of small-molecule BLI, the function is not intended to get accurate kinetics parameter, but to for the inital value to feed 'fitKinetics'
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitKinetics}}
}
\examples{
data.scaled=scaLoad(data.raw)
data.doubleRefed=doubleRef(data.scaled)
rate0=estKinetics(data.doubleRefed)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ dplot }% __ONLY ONE__ keyword per line
