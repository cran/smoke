\name{doubleRef}
\alias{doubleRef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
substract the double references
}
\description{
substract the double references from the BLI traces, align with the baselines before the association step; new dataframe returned as invisible and ploted
}
\usage{
doubleRef(data, tA0 = 1260, tD0 = 1860)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the BLI traces, first column is time in second
}
  \item{tA0}{
start-Time of Association step
}
  \item{tD0}{
start-Time of dissociation step
}
}
\details{
the input data may be raw data or the traces scaled with the protein-loading
}
\value{
dataframe: the first column is the time of association-then-dissociation segment; the next seven columns are the processed BLI traces
}
\author{
Qingan Sun, Xiaojun Li, James C Sacchettini
}
\examples{
data.scaled=scaLoad(data.raw)
data.doubleRefed=doubleRef(data.scaled)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
\keyword{ dplot }% __ONLY ONE__ keyword per line
