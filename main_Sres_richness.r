

library(nlme)
library(plot3D)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(directlabels)
library(MuMIn)
library(car)


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
m <- 1e-4
Kmean <- 10

###### Mono & Mixture: under frequent disturbance
T1 <- 50
T2 <- 0.5
T3 <- 50
dt <- 0.01


VfPert <- 0.
Valpha <- seq(0.2,0.8,length=7)
VrangeK <- seq(0.,0.7,length=8)
LL <- length(VfPert)*length(Valpha)*length(VrangeK)


###### parameters for scenario
    
noise.type <- 2     # 1 - random perturbations; 2 - predator perturbations
    
theta <- 1          # r-K trade-off: (-1)-positive correlation, (0)-independent; (1)-trade-off
gamma <- 0          # E-K trade-off: (-1)-positive correlation, (0)-independent; (1)-trade-off
    
    
nn_EFS <- c()

for(nn in 2:8)
{

  re_mono <- c()
  re_mix <- c()
  re_mix_sp <- c()
  
    for(k in 1:length(VfPert))
    {
      
      ###### Environmental (EE1) or predator perturbations (EE2)
      
      fPert <- VfPert[k]  ### fPert is the SD of noise (EE1) or the frequency of small perturbations (EE2)
    
      tt <- seq(0,T1,by=dt)
      
      eT1_1 <- matrix(rnorm(length(tt)*nn,0,fPert*4),nrow=nn)
      eT1_2 <- eT1_1*0
      eT1_2[,sample(1:length(tt),length(tt)*fPert)] <- 10
      
      eT2 <- matrix(rep(50,nn*T2/dt),nrow=nn)
      
      eT3 <- matrix(rep(0,nn*T3/dt),nrow=nn)
      
      
      for(kk in 1:length(VrangeK))
      {
        ###### species pool
        trait <- c(1:nn)
        VK <- seq(1-VrangeK[kk],1+VrangeK[kk],length=nn)*Kmean
        
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
        for(j in 1:nn)
        {
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
          resilience2 <- max(which(BT[-c(1:((T1+T2)/dt))] < (BT[(T1+T2)/dt]+(BT[length(BT)]-BT[(T1+T2)/dt])*0.5)))*dt

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
                                   BEf,BEs_temp,BEs_resist,BEs_resil,simpsonDiv,CE_SE,meanK))#,resistance.FI,resilience.FI,resilience.FI2))
          colnames(re_mix) <- c('noise','Alpha','deltaK','Richness','Functioning',
                                'Temporal stability','Resistance','Resilience','Species CV','Species Synchrony',
                                'BE_Func','BE_temp','BE_resist','BE_resil','Simpson Diversity',
                                'NBE','CE','SE','meanK')
        }
      }
    }
  nn_EFS <- rbind(nn_EFS,re_mix)
  print(nn)
}

#write.csv(nn_EFS, 'Sres_n=2_8.csv')



################################## Fig. 4a
##################################
nn2 <- nn_EFS
ind <- nn2[,4]%in%c(2:8) & nn2[,1]==0#& nn2[,2]%in%c(0.2,0.4,0.6,0.8) & nn2[,3]%in%c(0.1,0.3,0.5,0.7)

scatter3d(nn2[ind,5],nn2[ind,4],nn2[ind,7],group=as.factor(nn2[ind,4]),surface.alpha=0.2,
          axis.ticks=F,grid=T,grid.col='white',axis.scales = FALSE,surface.col=c(2,6,3,'orange',5,4,1),
          xlim=c(10,40),ylim=c(2,8),zlim=c(40,100),
          xlab='',ylab='',zlab='')

scatter3d(nn2[ind,5],nn2[ind,4],nn2[ind,8],group=as.factor(nn2[ind,4]),surface.alpha=0.2,
          axis.ticks=F,grid=T,grid.col='white',axis.scales = FALSE,surface.col=c(2,6,3,'orange',5,4,1),
          xlim=c(10,40),ylim=c(2,8),zlim=c(40,100),
          xlab='',ylab='',zlab='')


################################## Fig. 4b,c
##################################

### projection on funcitoning-stability space - different scenarios
scenEXT <- function(scenario1,nn_EFS){
  dd1 <- t(apply(scenario1,1,function(x){
    ind <- (abs(nn_EFS[,1]-x[1]) + abs(nn_EFS[,4]-x[2]) + abs(nn_EFS[,2]-x[3]) + abs(nn_EFS[,3]-x[4]))<1e-4;
    tmp <- if(!any(ind)){rep(NA,ncol(nn_EFS))}else{nn_EFS[ind,]}
    tmp}))
  dd1 <- dd1[!is.na(dd1[,1]),]
  return(dd1)
}

xylim <- function(zz){
  tmpx <- range(zz[,1])
  tmpy <- range(zz[,2])
  rx <- c((tmpx[1]*0.95),(tmpx[2]*1.05)); #rx[rx<0] <- rx[rx<0]/0.8*1.2
  ry <- c((tmpy[1]*0.95),(tmpy[2]*1.05)); #rx[rx<0] <- rx[rx<0]/0.8*1.2
  return(list(xx=rx,yy=ry))
}



par(mfcol=c(3,5),mar=c(3,3,1,1),mgp=c(1.5,0.5,0),tck=.02,lwd=1)
ccol <- c('red',rep('black',10))
ccol2 <- c('orange',rep('grey',10))
ccex <- c(2,rep(1,10))

###### variation in K increases with richness
nSp <- 2:8
ll <- length(nSp); aa <- 0.8

scenario1 <- cbind(noise.corr=0,richness=nSp,alpha=seq(aa,aa,length=ll),deltaK=seq(0.1,0.7,length=ll))
scenario2 <- cbind(noise.corr=0,richness=nSp,alpha=seq(aa,aa,length=ll),deltaK=seq(0.4,0.4,length=ll))
scenario3 <- cbind(noise.corr=0,richness=nSp,alpha=seq(aa,aa,length=ll),deltaK=seq(0.7,0.1,length=ll))

scenario1.2 <- cbind(noise.corr=0.5,richness=nSp,alpha=seq(aa,aa,length=ll),deltaK=seq(0.1,0.7,length=ll))
scenario2.2 <- cbind(noise.corr=0.5,richness=nSp,alpha=seq(aa,aa,length=ll),deltaK=seq(0.4,0.4,length=ll))
scenario3.2 <- cbind(noise.corr=0.5,richness=nSp,alpha=seq(aa,aa,length=ll),deltaK=seq(0.7,0.1,length=ll))

dd1 <- scenEXT(scenario1,nn_EFS)
dd2 <- scenEXT(scenario2,nn_EFS)
dd3 <- scenEXT(scenario3,nn_EFS)

dd1.2 <- scenEXT(scenario1.2,nn_EFS)
dd2.2 <- scenEXT(scenario2.2,nn_EFS)
dd3.2 <- scenEXT(scenario3.2,nn_EFS)

zz <- dd1[,c(5,7)]; zz0 <- zz
plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz0)$xx,ylim=xylim(zz0)$yy); lines(zz,col=1)
zz <- dd2[,c(5,7)]; zz0 <- zz
plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz0)$xx,ylim=xylim(zz0)$yy); lines(zz,col=1)
zz <- dd3[,c(5,7)]; zz0 <- zz
plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz0)$xx,ylim=xylim(zz0)$yy); lines(zz,col=1)


zz <- dd1[,c(5,8)]; zz0 <- zz
plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz0)$xx,ylim=xylim(zz0)$yy); lines(zz,col=1)
zz <- dd2[,c(5,8)]; zz0 <- zz
plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz0)$xx,ylim=xylim(zz0)$yy); lines(zz,col=1)
zz <- dd3[,c(5,8)]; zz0 <- zz
plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz0)$xx,ylim=xylim(zz0)$yy); lines(zz,col=1)


zz <- dd1[,17:18]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd2[,17:18]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd3[,17:18]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)

zz <- dd1[,c(4,19)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd2[,c(4,19)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd3[,c(4,19)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)

zz <- dd1[,c(4,15)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=c(0,8),ylim=c(0,8)); lines(zz,col=1); abline(0,1,lty=2)
zz <- dd2[,c(4,15)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=c(0,8),ylim=c(0,8)); lines(zz,col=1); abline(0,1,lty=2)
zz <- dd3[,c(4,15)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=c(0,8),ylim=c(0,8)); lines(zz,col=1); abline(0,1,lty=2)




###### competition increases with richness

par(mfcol=c(3,5),mar=c(2.5,2.5,1,1),mgp=c(1.2,0.2,0),tck=.02,lwd=1)
ccol <- c('red',rep('black',10))
ccex <- c(2,rep(1,10))

nSp <- 2:8;
ll<-length(nSp); dKK <- 0.5

scenario1 <- cbind(noise.corr=0,richness=nSp,alpha=seq(0.2,0.8,length=ll),deltaK=seq(dKK,dKK,length=ll))
scenario2 <- cbind(noise.corr=0,richness=nSp,alpha=seq(0.5,0.5,length=ll),deltaK=seq(dKK,dKK,length=ll))
scenario3 <- cbind(noise.corr=0,richness=nSp,alpha=seq(0.8,0.2,length=ll),deltaK=seq(dKK,dKK,length=ll))

scenario1.2 <- cbind(noise.corr=0.5,richness=nSp,alpha=seq(0.2,0.8,length=ll),deltaK=seq(dKK,dKK,length=ll))
scenario2.2 <- cbind(noise.corr=0.5,richness=nSp,alpha=seq(0.5,0.5,length=ll),deltaK=seq(dKK,dKK,length=ll))
scenario3.2 <- cbind(noise.corr=0.5,richness=nSp,alpha=seq(0.8,0.2,length=ll),deltaK=seq(dKK,dKK,length=ll))

dd1 <- scenEXT(scenario1,nn_EFS)
dd2 <- scenEXT(scenario2,nn_EFS)
dd3 <- scenEXT(scenario3,nn_EFS)

dd1.2 <- scenEXT(scenario1.2,nn_EFS)
dd2.2 <- scenEXT(scenario2.2,nn_EFS)
dd3.2 <- scenEXT(scenario3.2,nn_EFS)

zz <- dd1[,c(5,7)]; zz0 <- zz
plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz0)$xx,ylim=xylim(zz0)$yy); lines(zz,col=1)
zz <- dd2[,c(5,7)]; zz0 <- zz
plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz0)$xx,ylim=xylim(zz0)$yy); lines(zz,col=1)
zz <- dd3[,c(5,7)]; zz0 <- zz
plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz0)$xx,ylim=xylim(zz0)$yy); lines(zz,col=1)


zz <- dd1[,c(5,8)]; zz0 <- zz
plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz0)$xx,ylim=xylim(zz0)$yy); lines(zz,col=1)
zz <- dd2[,c(5,8)]; zz0 <- zz
plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz0)$xx,ylim=xylim(zz0)$yy); lines(zz,col=1)
zz <- dd3[,c(5,8)]; zz0 <- zz
plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz0)$xx,ylim=xylim(zz0)$yy); lines(zz,col=1)

zz <- dd1[,17:18]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd2[,17:18]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd3[,17:18]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)

zz <- dd1[,c(4,19)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd2[,c(4,19)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd3[,c(4,19)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)

zz <- dd1[,c(4,15)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=c(0,8),ylim=c(0,8)); lines(zz,col=1); abline(0,1,lty=2)
zz <- dd2[,c(4,15)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=c(0,8),ylim=c(0,8)); lines(zz,col=1); abline(0,1,lty=2)
zz <- dd3[,c(4,15)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=c(0,8),ylim=c(0,8)); lines(zz,col=1); abline(0,1,lty=2)

















##################################
##################################
ind <- nn_EFS[,4]%in%c(2,4,6,8)
scatter3d(nn_EFS[ind,5],nn_EFS[ind,4],nn_EFS[ind,7],group=as.factor(nn_EFS[ind,4]),surface.alpha=0.2,
          axis.ticks=T,grid=T,grid.col='grey',
          xlab='Functioning',ylab='Richness',zlab='Resistance')



##################################
##################################

### projection on funcitoning-stability space - different scenarios
scenEXT <- function(scenario1,nn_EFS){
  dd1 <- t(apply(scenario1,1,function(x){
    ind <- (abs(nn_EFS[,4]-x[1]) + abs(nn_EFS[,2]-x[2]) + abs(nn_EFS[,3]-x[3]))<1e-4;
    tmp <- if(!any(ind)){rep(NA,ncol(nn_EFS))}else{nn_EFS[ind,]}
    tmp}))
  dd1 <- dd1[!is.na(dd1[,1]),]
  return(dd1)
}

xylim <- function(zz){
  tmpx <- range(zz[,1])
  tmpy <- range(zz[,2])
  rx <- c((tmpx[1]*0.95),(tmpx[2]*1.05)); #rx[rx<0] <- rx[rx<0]/0.8*1.2
  ry <- c((tmpy[1]*0.95),(tmpy[2]*1.05)); #rx[rx<0] <- rx[rx<0]/0.8*1.2
  return(list(xx=rx,yy=ry))
}


##################################
###### variation in K increases with richness
##################################

par(mfcol=c(3,4),mar=c(3,3,1,1),mgp=c(1.5,0.5,0),tck=.02,lwd=1)
ccol <- c('red',rep('black',10))
ccex <- c(2,rep(1,10))

nSp <- 2:8
ll <- length(nSp); aa <- 0.4
scenario1 <- cbind(richness=nSp,alpha=seq(aa,aa,length=ll),deltaK=seq(0.1,0.7,length=ll))
scenario2 <- cbind(richness=nSp,alpha=seq(aa,aa,length=ll),deltaK=seq(0.4,0.4,length=ll))
scenario3 <- cbind(richness=nSp,alpha=seq(aa,aa,length=ll),deltaK=seq(0.7,0.1,length=ll))

dd1 <- scenEXT(scenario1,nn_EFS)
dd2 <- scenEXT(scenario2,nn_EFS)
dd3 <- scenEXT(scenario3,nn_EFS)

### resistance
zz <- dd1[,c(5,7)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd2[,c(5,7)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd3[,c(5,7)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)

### resilience
zz <- dd1[,c(5,8)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd2[,c(5,7)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd3[,c(5,7)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)

zz <- dd1[,17:18]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd2[,17:18]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd3[,17:18]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)

zz <- dd1[,c(4,15)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=c(0,8),ylim=c(0,8)); lines(zz,col=1); abline(0,1,lty=2)
zz <- dd2[,c(4,15)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=c(0,8),ylim=c(0,8)); lines(zz,col=1); abline(0,1,lty=2)
zz <- dd3[,c(4,15)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=c(0,8),ylim=c(0,8)); lines(zz,col=1); abline(0,1,lty=2)




##################################
###### competition increases with richness
##################################

par(mfcol=c(3,4),mar=c(2.5,2.5,1,1),mgp=c(1.2,0.2,0),tck=.02,lwd=1)
ccol <- c('red',rep('black',10))
ccex <- c(2,rep(1,10))

nSp <- 2:8;
ll<-length(nSp); dKK <- 0.1

scenario1 <- cbind(richness=nSp,alpha=seq(0.2,0.8,length=ll),deltaK=seq(dKK,dKK,length=ll))
scenario2 <- cbind(richness=nSp,alpha=seq(0.5,0.5,length=ll),deltaK=seq(dKK,dKK,length=ll))
scenario3 <- cbind(richness=nSp,alpha=seq(0.8,0.2,length=ll),deltaK=seq(dKK,dKK,length=ll))
dd1 <- scenEXT(scenario1,nn_EFS)
dd2 <- scenEXT(scenario2,nn_EFS)
dd3 <- scenEXT(scenario3,nn_EFS)

### resistance
zz <- dd1[,c(5,7)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd2[,c(5,7)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd3[,c(5,7)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)

### resilience
zz <- dd1[,c(5,8)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd2[,c(5,7)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd3[,c(5,7)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)

zz <- dd1[,17:18]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd2[,17:18]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)
zz <- dd3[,17:18]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=xylim(zz)$xx,ylim=xylim(zz)$yy); lines(zz,col=1)

zz <- dd1[,c(4,15)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=c(0,8),ylim=c(0,8)); lines(zz,col=1); abline(0,1,lty=2)
zz <- dd2[,c(4,15)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=c(0,8),ylim=c(0,8)); lines(zz,col=1); abline(0,1,lty=2)
zz <- dd3[,c(4,15)]; plot(zz,col=ccol,cex=ccex,pch=16,xlim=c(0,8),ylim=c(0,8)); lines(zz,col=1); abline(0,1,lty=2)







##################################

library("RColorBrewer")
colors <- brewer.pal(n=4, name="Dark2")

scatter3d(nn_EFS[,5],nn_EFS[,4],nn_EFS[,7],group=as.factor(nn_EFS[,4]),surface.alpha=0.2,
          surface.col=2:5,axis.ticks=T,grid=T,grid.col='grey',
          xlab='Functioning',ylab='Richness',zlab='Resistance')

scatter3d(nn_EFS[,5],nn_EFS[,4],nn_EFS[,8],group=as.factor(nn_EFS[,4]),
          surface.col=2:5,axis.ticks=T,grid=T,grid.col='white',
          xlab='Functioning',ylab='Richness',zlab='Resilience')

###### projection on funcitoning-stability space - different scenarios
scenEXT <- function(scenario,nn_EFS){
  dd1 <- t(apply(scenario1,1,function(x){
            ind <- nn_EFS[,4]==x[1] & nn_EFS[,2]==x[2] & nn_EFS[,3]==x[3];
            nn_EFS[ind,]}))
  return(dd1)
}

scenario1 <- cbind(richness=2:5,alpha=seq(0.5,0.5,length=4),deltaK=5)
dd1 <- scenEXT(scenario1,nn_EFS)
plot(dd1[,5],dd1[,7],type='b',xlab='Functioning',ylab='Resistance',pch=16,col=2:5)
plot(dd1[,5],dd1[,8],type='b',xlab='Functioning',ylab='Resilience',pch=16,col=2:5)

