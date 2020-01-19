setGeneric("estimate", function(obj) standardGeneric("estimate"))
setMethod("estimate", "Bli", function(obj) {
    t <- obj@traces[, 1]
    traces <- obj@traces[, 2:ncol(obj@traces)]
    tExp <- obj@tExp
    tAss <- tExp[2] - tExp[1]
    lig <- obj@lig
    if (length(lig) != ncol(traces)) 
        stop("the trace number mis-match ligand concentrations")
    totRange <- which(t >= tExp[1])
    tTot <- t[totRange]
    dTot <- traces[totRange, ]
    t <- tTot - tExp[1]
    tA <- t[which(t < tAss)]
    tD <- t[which(t > tAss)]
    dA <- dTot[which(t < tAss), ]
    dD <- dTot[which(t > tAss), ]
    rMax0 <- min(dTot)
    kOn <- kOff <- rep(NA, length(lig))
    for (i in 1:length(lig)) {
        y = dA[, i]
        aE = y[length(y)]
        yP = 1 - dA[, i]/aE
        tP = tA[which(yP > 0 & tA < (tAss/10))]
        yP = yP[which(yP > 0 & tA < (tAss/10))]
        yP = log(yP)
        modA = lm(yP ~ tP + 0)
        k0 = -coef(modA)
        modA = nls(y ~ A * (1 - exp(-k * tA)), start = list(k = k0, A = aE), trace = TRUE)
        kOn[i] = coef(modA)[1]
        y = dD[, i]
        d0 = predict(modA)[length(tA)]
        dE = y[length(y)]
        yP = (dD[, i] - dE)/(d0 - dE)
        tP = tD[which(yP > 0 & tD < (tAss * 1.1))] - 600
        yP = yP[which(yP > 0 & tD < (tAss * 1.1))]
        yP = log(yP)
        modD = lm(yP ~ tP + 0)
        k0 = -coef(modD)
        modD = nls(y ~ yE + (d0 - yE) * exp(-k * (tD - 600)), start = list(k = k0, 
            yE = dE))
        kOff[i] = coef(modD)[1]
    }
    obj@kOn0 <- mean(kOn/lig)
    obj@kOff0 <- mean(kOff)
    obj@status[4] <- TRUE
    return(obj)
})
