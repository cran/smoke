setGeneric("alignLoad", function(obj, loadStart, loadEnd) standardGeneric("alignLoad"))
setMethod("alignLoad", "Bli", function(obj, loadStart, loadEnd) {
    t <- obj@traces[, 1]
    traces <- obj@traces[, 2:ncol(obj@traces)]
    if (ncol(traces)%%2 != 0) 
        stop("the original traces should be paired")
    load0Range <- which(t < loadStart & t > (loadStart - 60))
    loadRange <- which(t < loadEnd & t > (loadEnd - 60))
    load0 <- load <- numeric()
    nSample <- ncol(traces)/2
    colnames(traces) <- 1:ncol(traces)
    for (i in 1:nSample) {
        tT <- t[load0Range]
        y <- traces[load0Range, i * 2 - 1]
        mod <- lm(y ~ tT)
        load0 <- c(load0, predict(mod, list(tT = loadStart)))
        y <- traces[load0Range, i * 2]
        mod <- lm(y ~ tT)
        load0 <- c(load0, predict(mod, list(tT = loadStart)))
        tT <- t[loadRange]
        y <- traces[loadRange, i * 2 - 1]
        mod <- lm(y ~ tT)
        load <- c(load, predict(mod, list(tT = loadEnd)))
        y <- traces[loadRange, i * 2]
        mod <- lm(y ~ tT)
        load <- c(load, predict(mod, list(tT = loadEnd)))
    }
    
    traces <- traces - rep(load0, each = length(t))
    load <- load - load0
    dim(load) <- c(2, nSample)
    meanLoad <- apply(load, 1, mean)
    loadFactor <- load/meanLoad
    loadFactor <- c(loadFactor)
    traces <- traces/rep(loadFactor, each = length(t))
    obj@traces <- data.frame(time = t, traces = traces)
    obj@status[1] <- TRUE
    return(obj)
})
