\name{data.raw}
\alias{data.raw}
\docType{data}
\title{
example raw data for BLI traces
}
\description{
example BLI traces read from csv file
}
\usage{data("data.raw")}
\format{
  A data frame with 12290 observations on the following 17 variables.
  \describe{
    \item{\code{t}}{a numeric vector}
    \item{\code{X1}}{a numeric vector}
    \item{\code{X2}}{a numeric vector}
    \item{\code{X3}}{a numeric vector}
    \item{\code{X4}}{a numeric vector}
    \item{\code{X5}}{a numeric vector}
    \item{\code{X6}}{a numeric vector}
    \item{\code{X7}}{a numeric vector}
    \item{\code{X8}}{a numeric vector}
    \item{\code{X9}}{a numeric vector}
    \item{\code{X10}}{a numeric vector}
    \item{\code{X11}}{a numeric vector}
    \item{\code{X12}}{a numeric vector}
    \item{\code{X13}}{a numeric vector}
    \item{\code{X14}}{a numeric vector}
    \item{\code{X15}}{a numeric vector}
    \item{\code{X16}}{a numeric vector}
  }
}
\details{
The BLI traces exported from Octet Red96 as SMPFormat, a specific text file. It is converted with python script into csv file first, then imported with readTraces function. The first column is time, the next 16 columns are time courses of BLI response.
}
\source{
Qingan Sun, Xiaojun Li, et al. The molecular basis of Pyrazinamide activity on Mycobacterium tuberculosis PanD. Nat Commun (in press)
}
\examples{
data(data.raw)
}
\keyword{datasets}
