setGeneric("fitKinetics", function(obj) standardGeneric("fitKinetics"))
setMethod("fitKinetics", "Bli", function(obj) {
    t <- obj@traces[, 1]
    traces <- obj@traces[, 2:ncol(obj@traces)]
    tExp <- obj@tExp
    tAss <- tExp[2] - tExp[1]
    lig <- obj@lig
    kOn0 <- obj@kOn0
    kOff0 <- obj@kOff0
    if ((length(kOn0) * length(kOff0)) == 0) 
        stop("no intial kOn0 and kOff0: input or estimate!")
    if (length(lig) != ncol(traces)) 
        stop("the trace number mis-match ligand concentrations")
    totRange <- which(t >= tExp[1])
    tTot <- t[totRange]
    dTot <- traces[totRange, ]
    t <- tTot - tExp[1]
    rMax0 <- min(dTot)
    x = rep(lig, each = length(tTot))
    t = rep(t, length(lig))
    y = unlist(dTot)
    temp = rep(0, length(lig))
    dLig = matrix(0, nrow = length(y), ncol = length(lig))
    iNum <- 1:length(lig)
    for (i in iNum) {
        temp1 = temp
        temp1[i] = 1
        dLig[, i] = rep(temp1, each = length(tTot))
    }
    mA = mD = rep(0, length(t))
    mA[which(t < tAss)] = 1
    mD[which(t > tAss)] = 1
    dfTot = data.frame(y, t, x, dLig, mA, mD)
    colnames(dfTot)[4:(3 + length(lig))] = paste("dLig", iNum, sep = "")
    startLst = c(rMax0, kOn0, kOff0, rep(0, 4 * length(lig)))
    names(startLst) <- c("rMax", "kOn", "kOff", paste("sa", iNum, sep = ""), paste("sd", 
        iNum, sep = ""), paste("ja", iNum, sep = ""), paste("jd", iNum, sep = ""))
    fmla <- "y ~ mA*(rMax/(1+(kOff/kOn/x))*(1-1/exp((kOn*x+kOff)*t))+"
    fmla <- paste(fmla, paste(paste(paste("dLig", iNum, sep = ""), paste(paste("ja", 
        iNum, sep = ""), paste(paste("sa", iNum, sep = ""), "t", sep = "*"), sep = "+"), 
        sep = "*("), collapse = ")+"), collapse = "")
    fmla <- paste(fmla, ")) + mD*( rMax/(1+(kOff/kOn/x))*(1-1/exp((kOn*x+kOff)*tAss))/exp(kOff*(t-tAss))+", 
        collapse = "")
    fmla <- paste(fmla, paste(paste(paste("dLig", iNum, sep = ""), paste(paste("jd", 
        iNum, sep = ""), paste(paste("sd", iNum, sep = ""), "(t-tAss)", sep = "*"), 
        sep = "+"), sep = "*("), collapse = ")+"), "))", collapse = "")
    fmla <- as.formula(fmla)
    modTot <- nls(fmla, data = dfTot, start = startLst, trace = TRUE)
    obj@kinetics <- modTot
    obj@status[5] <- TRUE
    return(obj)
})
