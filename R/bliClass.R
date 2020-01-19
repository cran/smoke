library(stats)
# Bli S4 class
setOldClass("nls")

Bli <- setClass("Bli", slots = c(traces = "data.frame", lig = "numeric", tExp = "numeric", 
    status = "vector", kinetics = "nls", kOn0 = "numeric", kOff0 = "numeric"))

setMethod("initialize", "Bli", function(.Object) {
    # .Object <- callNextMethod()
    .Object@status <- rep(FALSE, 5)
    names(.Object@status) <- c("alignLoad", "doubleBlank", "baseline", "estimate", 
        "fitKinetics")
    # .Object@kinetics <- structure(list(NULL), class = 'nls') .Object@kinetics <-
    # as(.Object@kinetics, 'nls') .Object <- callNextMethod()
    return(.Object)
})

setGeneric("traces", function(obj) standardGeneric("traces"))
setMethod("traces", "Bli", function(obj) {
    if (!length(obj@traces) > 0) 
        stop("please input traces in data.frame")
    return(obj@traces)
})

setGeneric("traces<-", function(obj, value) standardGeneric("traces<-"))
setMethod("traces<-", "Bli", function(obj, value) {
    if (!is.data.frame(value)) 
        stop("input traces in data.frame")
    obj@traces <- value
    return(obj)
})

setGeneric("ligand", function(obj) standardGeneric("ligand"))
setMethod("ligand", "Bli", function(obj) {
    if (!length(obj@lig) > 0) 
        stop("please input ligand concentration")
    return(obj@lig)
})

setGeneric("ligand<-", function(obj, value) standardGeneric("ligand<-"))
setMethod("ligand<-", "Bli", function(obj, value) {
    obj@lig <- value
    return(obj)
})

setGeneric("tExp", function(obj) standardGeneric("tExp"))
setMethod("tExp", "Bli", function(obj) {
    if (!length(obj@tExp) > 0) 
        stop("please input start Times of association and dissociation")
    return(obj@tExp)
})

setGeneric("tExp<-", function(obj, value) standardGeneric("tExp<-"))
setMethod("tExp<-", "Bli", function(obj, value) {
    obj@tExp <- value
    return(obj)
})

setGeneric("kOn0", function(obj) standardGeneric("kOn0"))
setMethod("kOn0", "Bli", function(obj) {
    return(obj@kOn0)
})

setGeneric("kOn0<-", function(obj, value) standardGeneric("kOn0<-"))
setMethod("kOn0<-", "Bli", function(obj, value) {
    obj@kOn0 <- value
    return(obj)
})

setGeneric("kOff0", function(obj) standardGeneric("kOff0"))
setMethod("kOff0", "Bli", function(obj) {
    return(obj@kOff0)
})

setGeneric("kOff0<-", function(obj, value) standardGeneric("kOff0<-"))
setMethod("kOff0<-", "Bli", function(obj, value) {
    obj@kOff0 <- value
    return(obj)
})

setGeneric("status", function(obj) standardGeneric("status"))
setMethod("status", "Bli", function(obj) {
    return(obj@status)
})

setGeneric("kinetics", function(obj) standardGeneric("kinetics"))
setMethod("kinetics", "Bli", function(obj) {
    if (!length(obj@kinetics) > 1) 
        stop("no fitting model yet")
    coefTot = coef(summary(obj@kinetics))
    kd = coefTot[3, 1]
    ka = coefTot[2, 1]
    kD = kd/ka
    kDsd = sqrt((coefTot[2, 2]/ka)^2 + (kd * coefTot[3, 2]/ka/ka)^2)
    coefTot = rbind(c(kD, kDsd, NA, NA), coefTot)
    rownames(coefTot)[1] = "KD"
    return(coefTot[1:4, ])
})
