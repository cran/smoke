\name{readTraces}
\alias{readTraces}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
read BLI traces from csv file
}
\description{
read BLI traces from csv file. The first column is time in second, the next 16 columns are BLI traces.The time series are put in an dataframe, returned as invisible and plot
}
\usage{
readTraces(csvFile = "traces.csv", tLoadS = 180, tLoadF = 780, tBlockS = 840, 
tBlockF = 960, tA0 = 1260, tD0 = 1860)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{csvFile}{
the csv file of BLI traces
}
  \item{tLoadS}{
start-Time of protein-Loading
}
  \item{tLoadF}{
end-Time of protein-Loading
}
  \item{tBlockS}{
start-Time of blocking
}
  \item{tBlockF}{
end-Time of blocking
}
  \item{tA0}{
start-Time of Association step
}
  \item{tD0}{
start-Time of dissociation step
}
}
\details{
The BLI traces exported from Octet Red96 as SMPFormat, a specific text file. It is converted with python script into csv file first. In this csv file, the first column is time, the next columns are BLI traces in the row-wise order.
}
\value{
dataframe: the first column is time; the next 16 columns are the raw BLI traces
}
\author{
Qingan Sun, Xiaojun Li, James C Sacchettini
}
\seealso{
\code{\link{data.raw}}
}
\examples{
\donttest{data.raw=readTraces("traces.csv")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
\keyword{ dplot }% __ONLY ONE__ keyword per line
