\name{fitKinetics}
\alias{fitKinetics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
non-linear regression of binding kinetics from BLI traces
}
\description{
fit the whole BLI dataset with one association-then-dissociation equation
}
\usage{
fitKinetics(data, tD0 = 600, kd0 = 0.08, ka0 = 0.013, lig = c(16/2^(0:6)), 
plotResidual = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the BLI traces with double-references corrected and aligned at start-time of association
}
  \item{tD0}{
start-Time of dissociation step
}
  \item{kd0}{
intial value of kOff, recommend to estimate with estKinetics function
}
  \item{ka0}{
intial value of kOn, recommend to estimate with estKinetics function
}
  \item{lig}{
concentrations of small-molecule analyte
}
  \item{plotResidual}{
logical, default as FALSE to plot the sensorGram with fitting; 'TRUE' to plot the residuals from fitting
}
}
\details{
To overcome the low SNR in the small-molecule BLI exepriment, the whole dataset is fitted with a single equation to minimize the model variability. The inital kOn and kOff values are recommended to estimate with 'estKinetics' function.
}
\value{
dataframe with the kinetics parameters and std.error from the non-linear regression,including KD, rMax, kOff and kOn
}
\author{
Qingan Sun, Xiaojun Li, James C Sacchettini
}
\note{
The initial kOn and kOff can be estimated with estKinetics function
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{estKinetics}}
}
\examples{
#data.raw=readTraces("traces.csv")
\donttest{data.scaled=scaLoad(data.raw)}
\donttest{data.doubleRefed=doubleRef(data.scaled)}
\donttest{rate0=estKinetics(data.doubleRefed)}
\donttest{kinetics=fitKinetics(data.doubleRefed)}
\donttest{kinetics=fitKinetics(data.doubleRefed,plotResidual=TRUE)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ dplot }% __ONLY ONE__ keyword per line
