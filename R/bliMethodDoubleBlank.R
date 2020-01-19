setGeneric("doubleBlank", function(obj) standardGeneric("doubleBlank"))
setMethod("doubleBlank", "Bli", function(obj) {
    if (obj@status[2]) 
        stop("already double-blank!")
    t <- obj@traces[, 1]
    traces <- obj@traces[, 2:ncol(obj@traces)]
    if (ncol(traces)%%2 != 0) 
        stop("the double-referenced traces should be paired")
    nSample <- ncol(traces)/2
    bkCor <- traces[, 2 * (1:nSample) - 1] - traces[, 2 * (1:nSample)]
    dbCor = bkCor[, 1:(nSample - 1)] - bkCor[, nSample]
    if (ncol(dbCor) != length(obj@lig)) 
        warning("ligand concentrations mis-match trace number") else colnames(dbCor) <- signif(obj@lig, 2)
    obj@traces <- data.frame(time = t, traces = dbCor)
    obj@status[2] <- TRUE
    return(obj)
})
