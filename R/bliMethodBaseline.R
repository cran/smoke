setGeneric("baseline", function(obj, tStart, tEnd) standardGeneric("baseline"))
setMethod("baseline", "Bli", function(obj, tStart, tEnd) {
    t <- obj@traces[, 1]
    traces <- obj@traces[, 2:ncol(obj@traces)]
    baseRange = which(t > tStart & t < tEnd)
    baseT = t[baseRange]
    baseline = numeric()
    for (i in 1:ncol(traces)) {
        y = traces[baseRange, i]
        mod = lm(y ~ baseT)
        baseline = c(baseline, predict(mod, list(baseT = tEnd)))
    }
    baseline = rep(baseline, each = length(t))
    dbNorm = traces - baseline
    if (ncol(dbNorm) != length(obj@lig)) 
        warning("ligand concentrations mis-match trace number") else colnames(dbNorm) <- signif(obj@lig, 2)
    obj@traces <- data.frame(time = t, traces = dbNorm)
    obj@status[3] <- TRUE
    return(obj)
})
