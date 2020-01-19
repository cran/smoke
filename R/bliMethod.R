setMethod("show", "Bli", function(object) {
    cat(is(object)[[1]], "Class \n==========\n")
    cat("status:\n")
    print(object@status)
    cat("--------------------\ntraces: data.frame with 1st column of time\n")
    if (length(object@traces) == 0) 
        cat("empty\n") else {
        print(object@traces[1:min(5, nrow(object@traces)), 1:min(6, ncol(object@traces))])
        cat(ncol(object@traces) - 1, "traces,", nrow(object@traces), "time points\n")
    }
    cat("--------------------\nligand: concentration from high to low\n")
    if (length(object@lig) == 0) 
        cat("empty\n") else print(object@lig)
    cat("--------------------\ntExp: start Time of associationa and dissociation\n")
    if (length(object@tExp) == 0) 
        cat("empty\n") else print(object@tExp)
    cat("--------------------\nkOn0: initial On-rate ")
    if (length(object@kOn0) == 0) 
        cat("(to be 'estimate'd or input)\n") else cat(object@kOn0, "\n")
    cat("--------------------\nkOn0: initial Off-rate ")
    if (length(object@kOff0) == 0) 
        cat("(to be 'estimate'd or input)\n") else cat(object@kOff0, "\n")
    cat("--------------------\nkinetics: fitting model\n")
    print(summary(object@kinetics))
    cat("--------------------\n")
})

setGeneric("plotTraces", function(obj, ...) standardGeneric("plotTraces"))
setMethod("plotTraces", "Bli", function(obj, ...) {
    if (!length(obj@traces) > 0) 
        stop("traces is empty")
    t <- obj@traces[, 1]
    traces <- obj@traces[, 2:ncol(obj@traces)]
    colFun = colorRampPalette(c("cyan", "magenta"))
    palet = colFun(ncol(traces))
    color = adjustcolor(palet, alpha.f = 0.3)
    plot(t, traces[, 1], type = "l", ylim = range(traces), col = color[1], xlab = "Time (s)", 
        ylab = "Response (nm)", ...)
    for (i in 2:ncol(traces)) lines(t, traces[, i], col = color[i], ...)
    abline(v = obj@tExp, lty = 2, ...)
})

setGeneric("plotKinetics", function(obj, ...) standardGeneric("plotKinetics"))
setMethod("plotKinetics", "Bli", function(obj, ...) {
    if (!length(obj@kinetics) > 1) 
        stop("not fitting kinetics yet")
    t <- obj@traces[, 1]
    traces <- obj@traces[, 2:ncol(obj@traces)]
    traces <- traces[which(t >= obj@tExp[1]), ]
    t <- t[which(t >= obj@tExp[1])] - obj@tExp[1]
    yF <- predict(obj@kinetics)
    dim(yF) <- dim(traces)
    colFun = colorRampPalette(c("cyan", "magenta"))
    color1 = colFun(ncol(traces))
    color = adjustcolor(color1, alpha.f = 0.3)
    plot(t, traces[, 1], type = "l", ylim = range(traces), col = color[1], xlab = "Time (s)", 
        ylab = "Response (nm)", ...)
    for (i in 2:ncol(traces)) lines(t, traces[, i], col = color[i], ...)
    for (i in 1:ncol(traces)) lines(t, yF[, i], col = color1[i], ...)
    abline(v = diff(obj@tExp), lty = 2, ...)
})

setGeneric("plotResiduals", function(obj, ...) standardGeneric("plotResiduals"))
setMethod("plotResiduals", "Bli", function(obj, ...) {
    if (!length(obj@kinetics) > 1) 
        stop("not fitting kinetics yet")
    t <- obj@traces[, 1]
    traces <- obj@traces[, 2:ncol(obj@traces)]
    traces <- traces[which(t >= obj@tExp[1]), ]
    t <- t[which(t >= obj@tExp[1])] - obj@tExp[1]
    yF <- predict(obj@kinetics)
    dim(yF) <- dim(traces)
    yF = traces - yF
    colFun = colorRampPalette(c("cyan", "magenta"))
    color1 = colFun(ncol(traces))
    color = adjustcolor(color1, alpha.f = 0.3)
    bound <- diff(range(traces))/2
    plot(t, yF[, 1], type = "l", ylim = c(-bound, bound), col = color[1], xlab = "Time (s)", 
        ylab = "Response (nm)", ...)
    for (i in 2:ncol(yF)) lines(t, yF[, i], col = color[i], ...)
    abline(v = diff(obj@tExp), lty = 2, ...)
    abline(h = 0, ...)
})

