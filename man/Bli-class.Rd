\name{Bli-class}
\Rdversion{1.1}
\docType{class}
\alias{Bli-class}

\title{Class \code{"Bli"}}
\description{
\code{Bli} is a S4 class to contain the Bli experiment data and analysis result.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{Bli(...)}.
}
\section{Slots}{
  \describe{
    \item{\code{traces}:}{Object of class \code{"data.frame"} of BLI traces; "time" in the 1st column }
    \item{\code{lig}:}{Object of class \code{"numeric"} of ligand concentrations }
    \item{\code{tExp}:}{Object of class \code{"numeric"} of association and dissociation times }
    \item{\code{status}:}{Object of class \code{"vector"} indicator and record of processing }
    \item{\code{kinetics}:}{Object of class \code{"nls"} model of binding kinetics }
    \item{\code{kOn0}:}{Object of class \code{"numeric"} initial value of on-rate }
    \item{\code{kOff0}:}{Object of class \code{"numeric"} initial value of off-rate }
  }
}
\section{Methods}{
  \describe{
    \item{alignLoad}{\code{signature(obj = "Bli")}: Align BLI traces with loading step }
    \item{baseline}{\code{signature(obj = "Bli")}: Align BLI traces with baseline step }
    \item{doubleBlank}{\code{signature(obj = "Bli")}: Substract the double-references }
    \item{estimate}{\code{signature(obj = "Bli")}: Estimate initial kinetic parameters }
    \item{fitKinetics}{\code{signature(obj = "Bli")}: Fit binding kinetics }
    \item{initialize}{\code{signature(.Object = "Bli")}: initializer }
    \item{kinetics}{\code{signature(obj = "Bli")}: Output kinetics model }
    \item{kOff0}{\code{signature(obj = "Bli")}: Getter for "Bli" slot of \code{kOff0} }
    \item{kOff0<-}{\code{signature(obj = "Bli")}: Setter for "Bli" slot of \code{kOff0} }
    \item{kOn0}{\code{signature(obj = "Bli")}: Getter for "Bli" slot of \code{kOn0} }
    \item{kOn0<-}{\code{signature(obj = "Bli")}: Setter for "Bli" slot of \code{kOn0} }
    \item{ligand}{\code{signature(obj = "Bli")}: Getter for \code{ligand} concentration }
    \item{ligand<-}{\code{signature(obj = "Bli")}: Setter for \code{ligand} concentration }
    \item{plotKinetics}{\code{signature(obj = "Bli")}: Plot kinetics model fitting }
    \item{plotResiduals}{\code{signature(obj = "Bli")}: Plot residuals from kinetics model fitting }
    \item{plotTraces}{\code{signature(obj = "Bli")}: Plot BLI traces }
    \item{show}{\code{signature(object = "Bli")}: show method }
    \item{status}{\code{signature(obj = "Bli")}: Getter for data processing \code{status} }
    \item{tExp}{\code{signature(obj = "Bli")}: Getter for BLI times \code{tExp} }
    \item{tExp<-}{\code{signature(obj = "Bli")}: Setter for BLI times \code{tExp} }
    \item{traces}{\code{signature(obj = "Bli")}: Getter of BLI \code{traces} }
    \item{traces<-}{\code{signature(obj = "Bli")}: Setter of BLI \code{traces} }
	 }
}
\references{
\insertRef{pza}{smoke}
}
\author{
\packageAuthor{smoke}

Maintainer: \packageMaintainer{smoke}
}
\examples{
showClass("Bli")
}
\keyword{classes}
\arguments{
  \item{\dots}{
.Object in \code{initialize}; leave empty
}
}
