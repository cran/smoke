#tLoadS=180;tLoadF=780   #Loading labelled Protein
#tBlockS=840;tBlockF=960 #Blocking with biocytin
#tA0=1260;tD0=tA0+600    #Starting of association and dissociation of Ligand

# read and plot the BLI traces from .csv file
readTraces=function(csvFile="traces.csv",tLoadS=180,tLoadF=780,tBlockS=840,tBlockF=960,tA0=1260,tD0=1860) {
   data=read.csv("traces.csv",header=FALSE)
   t=data[,1];raw=data[,1+(1:16)]
   colnames(raw)=1:16
   data=data.frame(t,raw)

   colFun=colorRampPalette(c("cyan","magenta"))
   palet=colFun(7)
   color1=c(palet,"black")
   color=adjustcolor(color1,alpha.f=0.3)
   plot(t,raw[,1],type='l',ylim=range(raw),col=color1[1],xlab="Time (s)",ylab="Response (nm)",lwd=8)
   lines(t,raw[,2],col=color1[1],lwd=8)
   for (i in 2:8) {lines(t,raw[,i*2-1],col=color1[i],lwd=8);lines(t,raw[,i*2],col=color1[i],lwd=8)}
   abline(v=c(tLoadS,tLoadF,tBlockS,tBlockF,tA0,tD0),lty=2,lwd=6)
 
   invisible(data)
}

#data.raw=readTraces()

# [option] scaling with Loading for samples and references separately
scaLoad=function(data,tLoadS=180,tLoadF=780) {
   t=data[,1];raw=data[,2:17]
   load0Range=which(t<tLoadS &t>(tLoadS-60));loadRange=which(t<tLoadF & t>(tLoadF-60))
   load0=load=numeric()
   for (i in 1:8) {
       tT=t[load0Range]
       y=raw[load0Range,i*2-1];mod=lm(y~tT);load0=c(load0,predict(mod,list(tT=tLoadS)))
       y=raw[load0Range,i*2];mod=lm(y~tT);load0=c(load0,predict(mod,list(tT=tLoadS)))
       tT=t[loadRange]
       y=raw[loadRange,i*2-1];mod=lm(y~tT);load=c(load,predict(mod,list(tT=tLoadF)))
       y=raw[loadRange,i*2];mod=lm(y~tT);load=c(load,predict(mod,list(tT=tLoadS)))
   }
   load0[(1:8)*2]=0;load[(1:8)*2]=1
   raw=raw-rep(load0,each=length(t))
   load=load-load0
   dim(load)=c(2,8);meanLoad=apply(load,1,mean)
   loadFactor=load/meanLoad;loadFactor=c(loadFactor)
   raw=raw/rep(loadFactor,each=length(t))

   colFun=colorRampPalette(c("cyan","magenta"))
   palet=colFun(7)
   color1=c(palet,"black")
   color=adjustcolor(color1,alpha.f=0.3)
   plot(t,raw[,1],type='l',ylim=range(raw),col=color1[1],xlab="Time (s)",ylab="Response (nm)",lwd=8)
   lines(t,raw[,2],col=color1[1],lwd=8)
   for (i in 2:8) {lines(t,raw[,i*2-1],col=color1[i],lwd=8);lines(t,raw[,i*2],col=color1[i],lwd=8)}
   abline(v=c(tLoadS,tLoadF),lty=2,lwd=6)

   data=data.frame(t,raw)
   invisible(data)
}

#data.scaled=scaLoad(data.raw)

doubleRef=function(data,tA0=1260,tD0=1860) {
   t=data[,1];raw=data[,2:17]
# substract the no-protein references
   bkCor = raw[,2*(1:8)-1]-raw[,2*(1:8)]
# substract the no-ligand reference
   dbCor=bkCor-bkCor[,8]

# align with the ending of baseline Step
   baseRange=which(t<tA0 & t>(tA0-180));baseT=t[baseRange]
   baseline=numeric()

   for (i in 1:8) {
     y=dbCor[baseRange,i]
     mod=lm(y~baseT)
     baseline=c(baseline,predict(mod,list(baseT=1260)))
   }
   baseline=rep(baseline,each=length(t))
   dbNorm=dbCor-baseline
   totRange=which(t>=tA0)
   t=t[totRange]-tA0;dTot=dbNorm[totRange,1:7]
   data=data.frame(t,dTot)

   colFun=colorRampPalette(c("cyan","magenta"))
   palet=colFun(7)
   color1=c(palet,"black")
   color=adjustcolor(color1,alpha.f=0.5)
   plot(t,dTot[,1],type='l',ylim=range(dTot),col=color1[1],
       xlab="Time (s)", ylab="Residuals (nm)")
   for (i in 2:7) lines(t,dTot[,i],col=color1[i])
   abline(v=c(tD0-tA0),lty=2)
   abline(h=0)
   invisible(data)
}

#data.doubleRefed=doubleRef(data.scaled)

# estimate kd0 and ka0 from data only (no start value)
estKinetics=function(data,tD0=600,lig=c(16/2^(0:6))) {
   t=data[,1];dTot=data[,2:8]
   tA=t[which(t<tD0)];tD=t[which(t>tD0)]
   dA=dTot[which(t<tD0),];dD=dTot[which(t>tD0),]
   rMax0=min(dTot);kOn=kOff=rep(NA,7)
   colFun=colorRampPalette(c("cyan","magenta"))
   palet=colFun(7)
   color1=c(palet,"black")
   color=adjustcolor(color1,alpha.f=0.5)
   plot(1,xlim=range(t),ylim=range(dTot),type='n',xlab="Time (s)", ylab="Response (nm)",col=color[1],lwd=6)
   for (i in 1:7) {
      y=dA[,i];aE=y[length(y)]
      lines(tA,dA[,i],col=color[i],lwd=6)
      yP=1-dA[,i]/aE
      tP=tA[which(yP>0 & tA<50)];yP=yP[which(yP>0 & tA<50)];yP=log(yP)
      modA=lm(yP~tP+0);k0=-coef(modA)
      modA=nls(y~A*(1-exp(-k*tA)),start=list(k=k0,A=aE),trace=TRUE)
      lines(tA,predict(modA),col=color1[i],lwd=6)
      kOn[i]=coef(modA)[1]
      y=dD[,i];d0=predict(modA)[length(tA)];dE=y[length(y)]
      lines(tD,dD[,i],col=color[i],lwd=6)
      yP=(dD[,i]-dE)/(d0-dE)
      tP=tD[which(yP>0 & tD<(tD0+50))]-tD0;yP=yP[which(yP>0 & tD <(tD0+50))];yP=log(yP)
      modD=lm(yP~tP+0);k0=-coef(modD)
      modD=nls(y~yE+(d0-yE)*exp(-k*(tD-tD0)),start=list(k=k0,yE=dE))
      lines(tD,predict(modD),col=color1[i],lwd=6)
      kOff[i]=coef(modD)[1]
   }
   abline(v=tD0,lwd=6,lty=2)

   kd0=mean(kOff); ka0=mean(kOn/lig)
   rate0=c(kd0,ka0); names(rate0)=c("kd0","ka0")
   invisible(rate0)
}

#rate0=estKinetics(data.doubleRefed)
#       kd0        ka0
# 0.08034611 0.01314539

fitKinetics=function(data,tD0=600,kd0=0.08,ka0=0.013,lig=c(16/2^(0:6)),plotResidual=FALSE) {
   t=tTot=data[,1];dTot=data[,2:8]
   rMax0=min(dTot)
   x=rep(lig,each=length(tTot));t=rep(t,7);y=unlist(dTot[,1:7])
   temp=rep(0,7);dt=matrix(0,nrow=length(y),ncol=7)
   for (i in 1:7) {
     temp1=temp;temp1[i]=1
     dt[,i]=rep(temp1,each=length(tTot))
   }
   mA=mD=rep(0,length(t));mA[which(t<tD0)]=1;mD[which(t>tD0)]=1
   dfTot=data.frame(y,t,x,dt,mA,mD)
   colnames(dfTot)[4:10]=c("dt1","dt2","dt3","dt4","dt5","dt6","dt7")
   startLst=list(rMax=rMax0,kd=kd0,ka=ka0,sa1=0,sa2=0,sa3=0,sa4=0,sa5=0,sa6=0,sa7=0,sd1=0,sd2=0,sd3=0,sd4=0,sd5=0,sd6=0,sd7=0,ja1=0,ja2=0,ja3=0,ja4=0,ja5=0,ja6=0,ja7=0,jd1=0,jd2=0,jd3=0,jd4=0,jd5=0,jd6=0,jd7=0)
   modTot=nls( y~mA*(rMax/(1+(kd/ka/x))*(1-1/exp((ka*x+kd)*t))+dt1*(ja1+sa1*t)+dt2*(ja2+sa2*t)+dt3*(ja3+sa3*t)+dt4*(ja4+sa4*t)+dt5*(ja5+sa5*t)+dt6*(ja6+sa6*t)+dt7*(ja7+sa7*t) )+mD*( rMax/(1+(kd/ka/x))*(1-1/exp((ka*x+kd)*tD0))/exp(kd*(t-tD0))+dt1*(jd1+sd1*(t-tD0))+dt2*(jd2+sd2*(t-tD0))+dt3*(jd3+sd3*(t-tD0))+dt4*(jd4+sd4*(t-tD0))+dt5*(jd5+sd5*(t-tD0))+dt6*(jd6+sd6*(t-tD0))+dt7*(jd7+sd7*(t-tD0)) ),
      data=dfTot,start=startLst,trace=TRUE)
   dfTotPre=cbind(dfTot,predict(modTot))
   coefTot=coef(summary(modTot))
   kd=coefTot[2,1];ka=coefTot[3,1]
   kD=kd/ka;kDsd=sqrt((coefTot[2,2]/ka)^2+(kd*coefTot[3,2]/ka/ka)^2)
   coefTot=rbind(c(kD,kDsd,NA,NA),coefTot)
   rownames(coefTot)[1]="KD";rownames(coefTot)[3:4]=c("kOff","kOn")
   colFun=colorRampPalette(c("cyan","magenta"))
   palet=colFun(7)
   color1=c(palet,"black")
   color=adjustcolor(color1,alpha.f=0.5)

   if (!plotResidual) {
   #png("sensorGram.png",width=3000,height=3000,pointsize=96)
   plot(tTot,dTot[,1],ylim=range(y),type='l',xlab="Time (s)", ylab="Response (nm)",col=color[1],lwd=6)
   for (i in 2:7) lines(tTot,dTot[,i],col=color[i],lwd=6)
   for (i in 1:7) lines(tTot,dfTotPre[which(dfTot[,3]==lig[i]),13],col=color1[i],lwd=6)
   abline(v=tD0,lwd=6,lty=2)    #dev.off()
   } else {
   residual=dfTotPre[,13]-dfTot[,1];bound=(range(y)[2]-range(y)[1])/2
   #png("residuals.png",width=3000,height=3000,pointsize=96)
   plot(tTot,residual[which(dfTot[,3]==lig[1])],ylim=c(-bound,bound),type='l',xlab="Time (s)", ylab=expression(Delta*"Response (nm)"),col=color[1],lwd=6)
for (i in 2:7) lines(tTot,residual[which(dfTot[,3]==lig[i])],col=color[i],lwd=6)
   abline(h=0,lwd=6)
   abline(v=tD0,lwd=6,lty=2)   #dev.off()
   }
   invisible(coefTot[1:4,])
}

#kinetics=fitKinetics(data.doubleRefed)
#kinetics=fitKinetics(data.doubleRefed,plotResidual=TRUE)


