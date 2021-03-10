

library(nlme)
library(plot3D)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(directlabels)
library(MuMIn)


###### functions
LH.partition <- function(mono,mix,ww=rep(1,length(mono)))
{
  SR <- length(mono)
  deltaY <- sum(mix) - sum(mono*ww)/sum(ww)
  
  RY.E <- ww/sum(ww)
  RY.O <- mix/mono
  
  deltaRY <- RY.O - RY.E
  
  CE <- SR*mean(deltaRY)*mean(mono)
  SE <- SR*cov(deltaRY,mono)*(SR-1)/SR
  
  c(deltaY,CE,SE)
}


###### parameters
g <- 0.01
m <- 1e-10
Kmean <- 10
nn <- 5

###### Mono & Mixture: under frequent disturbance
T1 <- 50
T2 <- 0.5
T3 <- 50
dt <- 0.01

VfPert <- 0. #c(0,0.1,0.2)  # 0.1
Valpha <- seq(0.2,0.8,length=7)
VrangeK <- seq(0,5,length=6)


###### parameters for scenario

noise.type <- 2     # 1 - random perturbations; 2 - predator perturbations

theta <- -1          # r-K trade-off: (-1)-positive correlation, (0)-independent; (1)-trade-off
gamma <- 0          # E-K trade-off: (-1)-positive correlation, (0)-independent; (1)-trade-off


re_mono <- c()
re_mix <- c()
re_mix_sp <- c()

for(k in 1:length(VfPert))
{
  
  ###### Environmental (EE1) or predator perturbations (EE2)
  
  fPert <- VfPert[k]  ### fPert is the SD of noise (EE1) or the frequency of small perturbations (EE2)

  tt <- seq(0,T1,by=dt)
  
  #eT1_1 <- matrix(rnorm(length(tt)*nn,0,fPert*4),nrow=nn)
  eT1_1 <- matrix(rnorm(length(tt)*nn,0,fPert*4),nrow=nn) + 1/3*matrix(1,nrow=nn,ncol=1)%*%rnorm(length(tt),0,fPert*4)
  #eT1_1 <- matrix(1,nrow=nn,ncol=1) %*% eT1_1[1,]
  eT1_2 <- eT1_1*0
  eT1_2[,sample(1:length(tt),length(tt)*fPert)] <- 10
  eT2 <- matrix(rep(50,nn*T2/dt),nrow=nn)
  eT3 <- matrix(rep(0,nn*T3/dt),nrow=nn)
  
  
  for(kk in 1:length(VrangeK))
  {
    ###### species pool
    trait <- c(1:nn)
    VK <- seq(Kmean-VrangeK[kk],Kmean+VrangeK[kk],length=nn)
    
    Vr <- 2*(10/VK)^theta
    VE <- 10*(VK/10)^gamma

    ### response to perturbation depends on species trait
    EE1 <- -cbind(eT1_1,eT2*matrix(rep(VE,ncol(eT2)),nrow=nn)*g,eT3)
    EE2 <- -cbind(eT1_2,eT2,eT3)*matrix(rep(VE,ncol(EE1)),nrow=nn)*g
    
    if(noise.type==1){
      EE <- EE1
    }else{
      EE <- EE2
    }
    
    ###### monoculture
    for(j in 1:nn){
      B.mono <- 1
      x.mono <- trait[j]
      r <- Vr[j]
      K <- VK[j]
      
      Bt <- B.mono
      for(i in 1:ncol(EE)){
        B1 <- Bt[length(Bt)]
        
        ff <- r*(1-B1/K)+EE[j,i]
        B2 <- pmax(B1 + (ff*B1 + m)*dt,0)
        
        Bt <- c(Bt,B2)
      }
      
      Bt.equ <- Bt[(T1/dt/4*3):(T1/dt)]
      meanBt_eqn <- mean(Bt.equ)
      cvBt_eqn <- sd(Bt.equ)/mean(Bt.equ)
      
      resis_minBt.value <- min(Bt[-(1:(T1/dt))])
      resistance2 <- Bt[(T1+T2)/dt]/meanBt_eqn
      #resis_minBt.when <- which.min(Bt[-(1:10001)])
      
      #resil_recover.time <- max(which(cumsum(Bt[-(1:10000)])/(1:(length(Bt)-10000)) < meanBt_eqn*0.5))
      #resil_recover.time <- max(which(Bt[-c(1:((T1+T2)/dt))] < Bt[length(Bt)]*0.6))
      resilience2 <- max(which(Bt[-c(1:((T1+T2)/dt))] < (Bt[(T1+T2)/dt]+(Bt[length(Bt)]-Bt[(T1+T2)/dt])*0.5)))*dt
      
      re_mono <- rbind(re_mono,c(fPert,NA,VrangeK[kk],x.mono,meanBt_eqn,cvBt_eqn,resistance2,resilience2))
    }
  
    ###### mixture: population dynamics
    for(kkk in 1:length(Valpha))
    {
      alpha <- Valpha[kkk]
      Malpha <- matrix(alpha,ncol=nn,nrow=nn); diag(Malpha) <- 1;
      
      B0 <- rep(1,nn)
      xx <- trait
      
      Bt <- t(as.matrix(B0))
      for(i in 1:ncol(EE))
      {
        B1 <- Bt[nrow(Bt),]
        
        ff <- Vr*(1-(B1%*%Malpha)/VK)+EE[,i]
        B2 <- pmax(B1 + (ff*B1 + m)*dt,0)
        
        Bt <- rbind(Bt,B2)
      }
      
      ###### [summary] mono vs. mixture
      BT <- rowSums(Bt)
      BT.equ <- BT[(T1/dt/4*3):(T1/dt)]; 
      meanBT_eqn <- mean(BT.equ)
      cvBT_eqn <- sd(BT.equ)/mean(BT.equ)
      
      Bt.equ <- Bt[(T1/dt/4*3):(T1/dt),]; covBt <- cov(Bt.equ)
      syn <- sum(covBt)/(sum(sqrt(diag(covBt))))^2
      cvSpecies <- cvBT_eqn/sqrt(syn)
      
      resis_minBT.value <- min(BT)
      resistance2 <- BT[(T1+T2)/dt]/meanBT_eqn
      resistance.FI <- meanBT_eqn/abs(meanBT_eqn-BT[(T1+T2)/dt])
      
      resil_recover.time <- max(which(BT[-c(1:((T1+T2)/dt))] < BT[length(BT)]*0.5))
      resilience2 <- max(which(BT[-c(1:((T1+T2)/dt))] < (BT[(T1+T2)/dt]+(BT[length(BT)]-BT[(T1+T2)/dt])*0.5)))*dt
      resilience.FI <- abs(meanBT_eqn-BT[(T1+T2)/dt])/abs(BT[(T1+T2+1)/dt]-BT[(T1+T2)/dt])
      resilience.FI2 <- abs(meanBT_eqn-BT[(T1+T2)/dt])/abs(BT[(T1+T2+2)/dt]-BT[(T1+T2)/dt])
      
      ### calculate Diversity, BEf, and BEs
      ind.mono <- re_mono[,1]==VfPert[k] & re_mono[,3]==VrangeK[kk]
      
      bio.mix <- colMeans(Bt[(T1/dt/2):(T1/dt),])
      simpsonDiv <- 1/sum((bio.mix/sum(bio.mix))^2)
      meanK <- sum(bio.mix*VK)/sum(bio.mix)
      
      bio.mono <- re_mono[ind.mono,5]
      CE_SE <- LH.partition(bio.mono,bio.mix)
      BEf <- meanBT_eqn/median(bio.mono)
      
      BEs_temp <- median(re_mono[ind.mono,6])/cvBT_eqn
      BEs_resist <- resistance2/median(re_mono[ind.mono,7])
      BEs_resil <- median(re_mono[ind.mono,8])/resilience2
      
      ### summary      
      re_mix_sp <- rbind(re_mix_sp,cbind(fPert,alpha,VrangeK[kk],1:nn,colMeans(Bt[(T1/dt/2):(T1/dt),])))
      re_mix <- rbind(re_mix,c(fPert,alpha,VrangeK[kk],nn,meanBT_eqn,1/cvBT_eqn,resistance2,1/resilience2,cvSpecies,syn,
                               BEf,BEs_temp,BEs_resist,BEs_resil,simpsonDiv,CE_SE,meanK,resistance.FI,resilience.FI,resilience.FI2))
    }
    print(c(k,kk,kkk))
  }
}


############################################################
###### Plot
############################################################

###### 
ind0 <- re_mix[,1]==0

PPI=1000 
INVERSE_TEXT_DIMENSION=8  
HEIGHT=PPI*INVERSE_TEXT_DIMENSION*0.5
WIDTH=HEIGHT*1.2 

p.names<-c("Functioning","CE","SE","Resistance","Resilience",'Diversity','meanK')
for (i in 1:7)
{
  kk <- c(5,17,18,7,8,15,19)[i]
  rr<-re_mix[ind0,kk]
  
  png(filename = paste('Fig_1.',p.names[i],'2.png'),width=WIDTH , height = HEIGHT , res = PPI)
  par(mgp=c(1.8,0.5,0))
  filled.contour(Valpha,VrangeK,matrix(re_mix[ind0,kk],nrow=length(Valpha)),
                 xlab=expression(alpha),ylab=expression(delta[K]),#main=p.names[i],
                 levels=seq(min(rr),max(rr),length=20),
                 col=colorRampPalette(c("#0000ff","#ffffff","#ff0000"))(21) ); 
  #lines(Valpha,(1-Valpha)/(4*Valpha+1)*10,lty=2,col='white')
  
  dev.off()
}



############################################################
###### Plot (OLD)
############################################################

###### Figure S2.
### diversity
par(mfrow=c(2,3),mgp=c(1.6,0.5,0),tck=-0.02)

contour(Valpha,VrangeK,matrix(re_mix[ind0,15],nrow=length(Valpha)),xlab='Alpha',ylab='delta_K'); #legend('topright','Diversity',bty='n',text.col='blue',cex=1.)
lines(Valpha,(1-Valpha)/(4*Valpha+1)*10,lty=2,col='blue')
contour(Valpha,VrangeK,matrix(re_mix[ind0,19],nrow=length(Valpha)),xlab='Alpha',ylab='delta_K'); #legend('topright','Diversity',bty='n',text.col='blue',cex=1.)
lines(Valpha,(1-Valpha)/(4*Valpha+1)*10,lty=2,col='blue')


###### Figure 1.
par(mfrow=c(2,3),mgp=c(1.6,0.5,0),tck=-0.02)
### functioning
contour(Valpha,VrangeK,matrix(re_mix[ind0,5],nrow=length(Valpha)),xlab='Alpha',ylab='delta_K'); #legend('topright','Functioning',bty='n',text.col='blue',cex=1.)
  lines(Valpha,(1-Valpha)/(4*Valpha+1)*10,lty=2,col='blue')
contour(Valpha,VrangeK,matrix(re_mix[ind0,17],nrow=length(Valpha)),xlab='Alpha',ylab='delta_K'); #legend('topright','CE',bty='n',text.col='blue',cex=1.)
  lines(Valpha,(1-Valpha)/(4*Valpha+1)*10,lty=2,col='blue')
contour(Valpha,VrangeK,matrix(re_mix[ind0,18],nrow=length(Valpha)),xlab='Alpha',ylab='delta_K'); #legend('topright','SE',bty='n',text.col='blue',cex=1.)
  lines(Valpha,(1-Valpha)/(4*Valpha+1)*10,lty=2,col='blue')

### Figure 1e-f. 
contour(Valpha,VrangeK,matrix(re_mix[ind0,7],nrow=length(Valpha)),xlab='Alpha',ylab='rangeK'); #legend('topright','Resistance',cex=1.)
contour(Valpha,VrangeK,matrix(re_mix[ind0,8],nrow=length(Valpha)),xlab='Alpha',ylab='rangeK'); #legend('topright','Resilience',cex=1.)



###### Figure 2b,c
### resistance
  bb <- coef(reg <- lm((zz<-re_mix[ind0,7])~(x1<-re_mix[ind0,17])+(x2<-re_mix[ind0,18]))); summary(reg)$r.
  zztmp <- outer(v1<-seq(min(x1),max(x1),length=10),
               v2<-seq(min(x2),max(x2),length=10),
               function(x1,x2){bb[1]+bb[2]*x1+bb[3]*x2})
  bb1 <- coef(reg <- lm((zz<-re_mix[ind0,7])~(x1<-re_mix[ind0,17])+(x2<-re_mix[ind0,18])+I(x2^2))); summary(reg)$r.
  zztmp1 <- outer(v1<-seq(min(x1),max(x1),length=10),
                  v2<-seq(min(x2),max(x2),length=10),
                  function(x1,x2){bb1[1]+bb1[2]*x1+bb1[3]*x2+bb1[4]*x2^2})
scatter3D(x1,x2,zz,col='black',pch=16,bty='g',cex=0.8,
          zlim=c(0.125,0.16),theta = 40, phi = 10, ticktype = "detailed",
          xlab='CE',ylab='SE',zlab='',main='Resistance',
          surf = list(x = v1, y = v2, z = zztmp1, col='black',
                      facets = NA))
### resilience
  bb <- coef(reg <- lm((zz<-re_mix[ind0,8])~(x1<-re_mix[ind0,17])+(x2<-re_mix[ind0,18]))); summary(reg)$r.
  zztmp <- outer(v1<-seq(min(x1),max(x1),length=10),
               v2<-seq(min(x2),max(x2),length=10),
               function(x1,x2){bb[1]+bb[2]*x1+bb[3]*x2})
  bb1 <- coef(reg <- lm((zz<-re_mix[ind0,8])~(x1<-re_mix[ind0,17])+(x2<-re_mix[ind0,18])+I(x2^2))); summary(reg)$r.
  zztmp1 <- outer(v1<-seq(min(x1),max(x1),length=10),
                  v2<-seq(min(x2),max(x2),length=10),
                  function(x1,x2){bb1[1]+bb1[2]*x1+bb1[3]*x2+bb1[4]*x2^2})



### Figure 3c-f. 

###### statistics II: resistance vs. functioning 

par(mfrow=c(2,3),mgp=c(1.6,0.5,0),tck=-0.02)
# raw data
plot(re_mix[ind0,7]~re_mix[ind0,5],pch=16,log='',col='darkgrey',
     xlab='Functioning',ylab='Resistance')

# rangeK as a random factor
plot(re_mix[ind0,7]~re_mix[ind0,5],pch=16,log='',col=as.factor(re_mix[ind0,3]),
     xlab='Functioning',ylab='Resistance')

library(nlme)
rr <- lme(V7~V5,~1|V3,data=as.data.frame(re_mix[ind0,]))
abline(summary(rr)$coef$fixed,lwd=3)
for(k in 1:length(VrangeK)){
  ab <- summary(rr)$coef$fixed
  ab[1] <- ab[1]+summary(rr)$coef$random$V3[k]
  abline(ab,col=k,lty=2)
}

# alpha as a random factor
plot(re_mix[ind0,7]~re_mix[ind0,5],pch=16,log='',col=as.factor(re_mix[ind0,2]),
     xlab='Functioning',ylab='Resistance')

library(nlme)
rr <- lme(V7~V5,~1|V2,data=as.data.frame(re_mix[ind0,]))
abline(summary(rr)$coef$fixed,lwd=3)
for(k in 1:length(Valpha)){
  ab <- summary(rr)$coef$fixed
  ab[1] <- ab[1]+summary(rr)$coef$random$V2[k]
  abline(ab,col=k,lty=2)
}


###### statistics III: resilience vs. functioning 

# raw data
plot(re_mix[ind0,8]~re_mix[ind0,5],pch=16,log='',col='darkgrey',
     xlab='Functioning',ylab='Resilience')

# rangeK as a random factor
plot(re_mix[ind0,8]~re_mix[ind0,5],pch=16,log='',col=as.factor(re_mix[ind0,3]),
     xlab='Functioning',ylab='Resilience')

library(nlme)
rr <- lme(V8~V5,~1|V3,data=as.data.frame(re_mix[ind0,]))
abline(summary(rr)$coef$fixed,lwd=3)
for(k in 1:length(VrangeK)){
  ab <- summary(rr)$coef$fixed
  ab[1] <- ab[1]+summary(rr)$coef$random$V3[k]
  abline(ab,col=k,lty=2)
}

# alpha as a random factor
plot(re_mix[ind0,8]~re_mix[ind0,5],pch=16,log='',col=as.factor(re_mix[ind0,2]),
     xlab='Functioning',ylab='Resilience')

library(nlme)
rr <- lme(V8~V5,~1|V2,data=as.data.frame(re_mix[ind0,]))
abline(summary(rr)$coef$fixed,lwd=3)
for(k in 1:length(Valpha)){
  ab <- summary(rr)$coef$fixed
  ab[1] <- ab[1]+summary(rr)$coef$random$V2[k]
  abline(ab,col=k,lty=2)
}


################################################################
###### check metrics of resilience and resistance
################################################################

Bt1 <- Bt   ### a simulation with delta_K=0
Bt2 <- Bt   ### a simulation with delta_K=5
  
par(mfcol=c(1,1),mar=c(3,3,1,1),mgp=c(1.6,0.5,0))
plot.ts(EE[1,],ylab='EE',xlim=c(0,50),col='blue')
plot.ts(BT <- rowSums(Bt1),ylab='Biomass',log='y',xlim=c(4900,5400),ylim=c(0.2,20),lwd=2,col='red')
matplot(Bt1,type='l',log='',add=T,lwd=0.5)
abline(h=sum(Bt1[4900,]*0.01+Bt1[5052,]*0.99),col='red')

lines(rowSums(Bt2),lwd=2,lty=1)
matplot(Bt2,type='l',log='',add=T,lwd=2,lty=2,col='black')
abline(h=sum(Bt2[4900,]*0.01+Bt2[5052,]*0.99),col='black')

aaa <- Bt1[-c(1:5051),];
matplot((aaa[-1,]-aaa[-nrow(aaa),]),type='l',xlim=c(0,600),ylim=c(0,0.1))
lines(rowSums(aaa[-1,]-aaa[-nrow(aaa),]),col='red',lwd=2)
bbb <- Bt2[-c(1:5051),]
matplot((bbb[-1,]-bbb[-nrow(bbb),]),type='l',xlim=c(0,600),add=T,col='grey')
lines(rowSums(bbb[-1,]-bbb[-nrow(bbb),]),col=1,lwd=2)

### rescaled (see Ingrisch1 & Bahn 2018)
aaa2 <- rowSums(aaa);aaa2 <- (aaa2-min(aaa2))/(max(aaa2)-min(aaa2))
plot.ts(aaa2,xlim=c(0,500),col='red');abline(h=1,lty=2)
bbb2 <- rowSums(bbb);bbb2 <- (bbb2-min(bbb2))/(max(bbb2)-min(bbb2))
lines(bbb2,xlim=c(0,500),col=1)

