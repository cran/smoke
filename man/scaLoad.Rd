\name{scaLoad}
\alias{scaLoad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scaling BLI traces with protein-Loadings
}
\description{
An optional step in data processing
}
\usage{
scaLoad(data, tLoadS = 180, tLoadF = 780)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
the BLI traces, first column is time in second
}
  \item{tLoadS}{
start-Time of protein-Loading step
}
  \item{tLoadF}{
end-Time of protein-Loading step
}
}
\details{
This is an optional step in data processing. The BLI traces are scaled with protein-Loadings, the new dataframe returned as invisible and ploted.
}
\value{
dataframe: the first column is time; the next 16 columns are BLI traces
}
\author{
Qingan Sun, Xiaojun Li, James C Sacchettini
}
\examples{
data.scaled=scaLoad(data.raw)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
\keyword{ dplot }% __ONLY ONE__ keyword per line
