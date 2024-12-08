library(latex2exp)

tess<-function(V,S,K){#function to calculate tessarine coefficients
  eq1<-c(0,0,0,0,S^2,0,45*V^2-9*K,0,-36*V,0,36)
  if(V==0){
    return(c(0,0,0,1))
  }
  sols<-polyroot(eq1)
  ds<- Re(sols[which(abs(Im(sols))<10^(-4) & Re(sols)!=0)])
  if(length(ds)==0){
        fn<-function(theta,V,S,K){
        b<-theta[1]
        c<-theta[2]
        d<-theta[3]
        f1<-(V-b^2+c^2-d^2)
        f2<-(S-6*b*c*d)/3
        f3<-(K+b^4+c^4+d^4-6*V*(b^2-c^2+d^2)-6*(b^2*c^2-b^2*d^2+c^2*d^2))/12
        return((f1^2+f2^2+f3^2))
      }
      gr<-function(theta,V,S,K){
        b<-theta[1]
        c<-theta[2]
        d<-theta[3]
        fb<--4/3*c*d*(S-6*b*c*d)+4*b*(b^2-c^2+d^2-V)+1/18*b*(b^2-3*(c^2-d^2+V))*(b^4+c^4+d^4+K-6*d^2*V+6*c^2*(V-d^2)-6*b^2*(c^2-d^2+V))
        fc<--4/3*b*d*(S-6*b*c*d)+4*c*(V-b^2+c^2-d^2)+1/18*c*(-3*b^2+c^2-3*d^2+3*V)*(b^4+c^4+d^4+K-6*d^2*V+6*c^2*(V-d^2)-6*b^2*(c^2-d^2+V))
        fd<--4/3*b*c*(S-6*b*c*d)+4*d*(b^2-c^2+d^2-V)+1/18*d*(3*b^2-3*c^2+d^2-3*V)*(b^4+c^4+d^4+K-6*d^2*V+6*c^2*(V-d^2)-6*b^2*(c^2-d^2+V))
        return(c(fb,fc,fd))
      }
      temp<-function(sdw){
        V<-V/sdw^2
        S<-S/sdw^3
        K<-K/sdw^4
        f1<-function(b,c,d){
          return(V-b^2+c^2-d^2)
        }
        f2<-function(b,c,d){
          return(S-6*b*c*d)
        }
        f3<-function(b,c,d){
          return(K+b^4+c^4+d^4-6*V*(b^2-c^2+d^2)-6*(b^2*c^2-b^2*d^2+c^2*d^2))
        }
        out<-optim(c(sqrt(V),.1,-.1),fn=fn,gr=gr,method="BFGS",V=V,S=S,K=K)$par
        
        b<-out[1]
        c<-out[2]
        d<-out[3]
        err<-abs(f1(b,c,d))/V+abs(f2(b,c,d)/3/S)+abs(f3(b,c,d)/6)/K
        return(err)
      }
      sdw<-optim(1,fn=temp,method="Brent",lower=.3,upper=5)$par
      out<-optim(c(sqrt(V),.1,-.1),fn=fn,gr=gr,method="BFGS",V=V/sdw^2,S=S/sdw^3,K=K/sdw^4)$par
    }
  if(length(ds)!=0){
  min<-1000000000000000
  for(d in ds){
    eq2<-c(-S^2,0,36*d^2*(V-d^2),0,36*d^2)
    sols2<-polyroot(eq2)
    cs<- Re(sols2[which(abs(Im(sols2))<10^(-4))])
    if(length(cs)==0){
      break
    }
    for(c in cs){
      b<-sign(S)/sign(c)/sign(d)*sqrt(V+c^2-d^2)
      temp<-sqrt(b^2+c^2+d^2)
      if(temp< min){
        sdw<-1
        out<-c(b,c,d)
        min<-temp
      }
    }
  }
  }
  if(length(cs)==0||length(ds)==0){
    fn<-function(theta,V,S,K){
      b<-theta[1]
      c<-theta[2]
      d<-theta[3]
      f1<-(V-b^2+c^2-d^2)
      f2<-(S-6*b*c*d)/3
      f3<-(K+b^4+c^4+d^4-6*V*(b^2-c^2+d^2)-6*(b^2*c^2-b^2*d^2+c^2*d^2))/12
      return((f1^2+f2^2+f3^2))
    }
    gr<-function(theta,V,S,K){
      b<-theta[1]
      c<-theta[2]
      d<-theta[3]
      fb<--4/3*c*d*(S-6*b*c*d)+4*b*(b^2-c^2+d^2-V)+1/18*b*(b^2-3*(c^2-d^2+V))*(b^4+c^4+d^4+K-6*d^2*V+6*c^2*(V-d^2)-6*b^2*(c^2-d^2+V))
      fc<--4/3*b*d*(S-6*b*c*d)+4*c*(V-b^2+c^2-d^2)+1/18*c*(-3*b^2+c^2-3*d^2+3*V)*(b^4+c^4+d^4+K-6*d^2*V+6*c^2*(V-d^2)-6*b^2*(c^2-d^2+V))
      fd<--4/3*b*c*(S-6*b*c*d)+4*d*(b^2-c^2+d^2-V)+1/18*d*(3*b^2-3*c^2+d^2-3*V)*(b^4+c^4+d^4+K-6*d^2*V+6*c^2*(V-d^2)-6*b^2*(c^2-d^2+V))
      return(c(fb,fc,fd))
    }
    temp<-function(sdw){
      V<-V/sdw^2
      S<-S/sdw^3
      K<-K/sdw^4
      f1<-function(b,c,d){
        return(V-b^2+c^2-d^2)
      }
      f2<-function(b,c,d){
        return(S-6*b*c*d)
      }
      f3<-function(b,c,d){
        return(K+b^4+c^4+d^4-6*V*(b^2-c^2+d^2)-6*(b^2*c^2-b^2*d^2+c^2*d^2))
      }
      out<-optim(c(sqrt(V),.1,-.1),fn=fn,gr=gr,method="BFGS",V=V,S=S,K=K)$par
      b<-out[1]
      c<-out[2]
      d<-out[3]
      if(Skew==0){
        err<-abs(f1(b,c,d)/V+abs(f3(b,c,d)/6/K))
      }else{
        err<-abs(f1(b,c,d))/V+abs(f2(b,c,d)/3/S)+abs(f3(b,c,d)/6)/K
      }
      return(err)
    }
    sdw<-optim(1,fn=temp,method="Brent",lower=.3,upper=5)$par
    out<-optim(c(sqrt(V),.1,-.1),fn=fn,gr=gr,method="BFGS",V=V/sdw^2,S=S/sdw^3,K=K/sdw^4)$par
  }
  return(c(out,sdw))
}
err<-function(b,c,d){#remainder of bias reduction
  f1<-function(b,c,d){
    return(Va/sdw^2-b^2+c^2-d^2)
  }
  f2<-function(b,c,d){
    return(Skew/sdw^3-6*b*c*d)
  }
  f3<-function(b,c,d){
    return(Kurt/sdw^4+b^4+c^4+d^4-6*Va*(b^2-c^2+d^2)-6*(b^2*c^2-b^2*d^2+c^2*d^2))
  }
  if(Skew==0){
    err<-abs(f1(b,c,d)/Va*sdw^2+abs(f3(b,c,d)/6/Kurt*sdw^4))
  }else{
    err<-abs(f1(b,c,d))/Va*sdw^2+abs(f2(b,c,d)/3/Skew*sdw^3)+abs(f3(b,c,d)/6)/Kurt*sdw^4
  }
  return(err)
}

sufflog<-function(data){#Sufficiency estimator
  w<-data[1:n]-MU
  y<-data[(n+1):(2*n)]
  muw<-0#mean(w)
  sdw<-1#sqrt(var(w))
  score<-function(theta){
    b0<-theta[1]
    b1<-theta[2]
    deltastar<-w+y*Va*b1-1/2*Va*b1
    score0<-sum(y-1/(1+exp(-b0-b1*deltastar)))^2
    score1<-sum((y-1/(1+exp(-b0-b1*deltastar)))*deltastar)^2
    return(score0+score1)
  }
  grad<-function(theta){
    b0<-theta[1]
    b1<-theta[2]
    deltastar<-w+y*Va*b1-1/2*Va*b1
    df1b0<-1/(-2-2*cosh(b0+b1*w+1/2*b1^2*Va*(2*y-1)))
    df1b1<- -(w+b1*Va*(2*y-1))/(2*(1+cosh(b0+b1*w+1/2*b1^2*Va*(2*y-1))))
    df2b0<- (b1*Va-2*w-2*b1*Va*y)/(4+4*cosh(b0+b1*w+1/2*b1^2*Va*(2*y-1)))
    df2b1<-Va*(y-1/2)*(y-1+1/(1+exp(b0+b1*(w+b1*Va*(y-1/2)))))-1/8*(w+b1*Va*(2*y-1))*(2*w+b1*Va*(2*y-1))*1/cosh(1/4*(2*b0+b1*(2*w+b1*Va*(2*y-1))))^2
    f1<-(y-1/(1+exp(-b0-b1*deltastar)))
    f2<-((y-1/(1+exp(-b0-b1*deltastar)))*deltastar)
    gr1<-2*mean(f1)*mean(df1b0)+2*mean(f2)*mean(df2b0)
    gr2<-2*mean(f1)*mean(df1b1)+2*mean(f2)*mean(df2b1)
    return(c(gr1,gr2))
  }
  init<-naivelog(data)
  beta<-optim(init,fn=score,method="BFGS",gr=grad)$par
  b1<-beta[2]
  b0<-beta[1]
  return(c(b0,b1))
}
mccs<-function(data){#Monte-Carlo Corrected score
  w<-data[1:n]-MU
  y<-data[(n+1):(2*n)]
  muw<-0#mean(w)
  sdw<-1#sqrt(var(w))
  wstar<-(w-muw)/sdw
  Vastar<-Va/sdw^2
  score<-function(theta){
    b0<-theta[1]
    b1<-theta[2]
    Splus1<-0
    Sminus1<-0
    Splus2<-0
    Sminus2<-0
    K<-1:20
    for(k in K){
      Splus1<-Splus1+(-1)^(k+1)*exp(-k*(b0+b1*wstar)-k^2*b1^2/2*Vastar)
      Sminus1<-Sminus1+(-1)^(k+1)*exp(k*(b0+b1*wstar)-k^2*b1^2/2*Vastar)
      Splus2<-Splus2+(-1)^(k+1)*exp(-k*(b0+b1*wstar)-k^2*b1^2/2*Vastar)*(wstar+k*b1*Vastar)
      Sminus2<-Sminus2+(-1)^(k+1)*exp(k*(b0+b1*wstar)-k^2*b1^2/2*Vastar)*(wstar-k*b1*Vastar)
    }
    A1<-is.nan(I(b0+b1*wstar<0)*Sminus1)
    A2<-is.nan(I(b0+b1*wstar>0)*(1-Splus1))
    B1<-is.nan(I(b0+b1*wstar<0)*Sminus2)
    B2<-is.nan(I(b0+b1*wstar>0)*(wstar-Splus2))
    Sminus1[A1]<-0
    Splus1[A2]<-0
    Sminus2[B1]<-0
    Splus2[B2]<-0
    score1<-mean(y-I(b0+b1*wstar<0)*Sminus1-I(b0+b1*wstar>0)*(1-Splus1))
    score2<-mean(y*wstar-I(b0+b1*wstar<0)*Sminus2-I(b0+b1*wstar>0)*(wstar-Splus2))
    return((score1)^2+(score2)^2)
  }
  init<-naivelog(data)
  beta1<-optim(init,fn=score,method="BFGS")$par
  beta2<-optim(init,fn=score,method="Nelder-Mead")$par
  if(score(beta1)<score(beta2)){
    beta<-beta1
  }else{
    beta<-beta2
  }
  b1<-beta[2]
  b0<-beta[1]
  return(beta)
}


naivelog<-function(data){#Naive estimator
  w<-data[1:n]-MU
  y<-data[(n+1):(2*n)]
  muw<-0
  sdw<-1#sqrt(var(w))
  wstar<-(w-muw)/sdw
  score<-function(theta){
    b0<-theta[1]
    b1<-theta[2]
    score1<-(mean(y*wstar-wstar/(1+exp(-b0-b1*wstar))))^2
    score0<-(mean(y-1/(1+exp(-b0-b1*wstar))))^2
    return(score0+score1)
  }
  Hess<-function(theta){
    b0<-theta[1]
    b1<-theta[2]
    phi1b0<- -2*(sum(y*wstar)-sum(wstar/(1+exp(-b0-b1*wstar))))*sum(wstar/(2+2*cosh(b0+b1*wstar)))/n^2
    phi1b1<- -2*(sum(y*wstar)-sum(wstar/(1+exp(-b0-b1*wstar))))*sum(exp(-b0-b1*wstar)*wstar^2/(1+exp(-b0-b1*wstar))^2)/n^2
    phi2b0<--2*(sum(y)-sum(1/(1+exp(-b0-b1*wstar))))*sum(1/(2+2*cosh(b0+b1*wstar)))/n^2
    phi2b1<- -2*(sum(y)-sum(1/(1+exp(-b0-b1*wstar))))*sum(wstar/(2+2*cosh(b0+b1*wstar)))/n^2
    return(c(phi1b0+phi2b0,phi1b1+phi2b1))
  }
  beta1<-optim(c(0,0),fn=score,gr=Hess,method="Nelder-Mead")$par
  return(beta1)
}

newlog<-function(data){#novel method
  w<-data[1:n]-MU
  y<-data[(n+1):(2*n)]
  muw<-0#mean(w)
  wstar<-(w-muw)/sdw
  sdw<-1#sqrt(var(w))
  score<-function(theta){
    b0<-theta[1]
    b1<-theta[2]
    score0<-mean((y-Re((1+1/2*exp(-b0-b1*(b*1i+c+d*1i+wstar))*(1+exp(2*b1*(c+d*1i))))/(1+exp(-2*b0+2*b1*(c+d*1i)-2*b1*(b*1i+c+d*1i+wstar))+exp(-b0-b1*(b*1i+c+d*1i+wstar))+exp(-b0+2*b1*(c+d*1i)-b1*(b*1i+c+d*1i+wstar))))))^2
    score1<-mean(y*wstar-Re(exp(b0+b1*wstar+b1*b*1i)*(exp(b0+b1*wstar+b1*b*1i)*(wstar+b*1i)+(wstar+b*1i)*cosh(b1*(c+d*1i))+(c+d*1i)*sinh(b1*(c+d*1i)))/(1+exp(2*(b0+b1*b*1i+b1*wstar))+2*exp(b0+b1*wstar+b1*b*1i)*cosh(b1*(c+d*1i)))))^2
    return(score0+score1)
  }
  init<-naivelog(data)
  beta1<-optim(init,fn=score,method="BFGS")$par
  beta2<-optim(init,fn=score,method="Nelder-Mead")$par
  if(score(beta1)<score(beta2)){
    beta<-beta1
  }else{
    beta<-beta2
  }
  return(beta)
}
set.seed(1)
reliability<-seq(.8,1,.1)
samp<-1000
n<-1000
set.seed(1)
newbias1<-c()
newmse1<-c()
naibias1<-c()
naimse1<-c()
commse1<-c()
combias1<-c()
newbias2<-c()
newmse2<-c()
naibias2<-c()
naimse2<-c()
commse2<-c()
combias2<-c()
mccsbias1<-c()
mccsmse1<-c()
mccsbias2<-c()
mccsmse2<-c()
suffbias1<-c()
suffmse1<-c()
suffbias2<-c()
suffmse2<-c()
newiqr1<-c()
newiqr2<-c()
suffiqr1<-c()
suffiqr2<-c()
naiiqr1<-c()
naiiqr2<-c()
mccsiqr1<-c()
mccsiqr2<-c()
x<-matrix(rnorm(n*samp,1,1),nrow=n,ncol=samp)
b0<-5
b1<--5
pk<-1/(1+exp(-b0-b1*x))
y<-apply(pk,2,rbinom,n=n,size=1)
MU<-0
V<-var(matrix(x,ncol=1))[1]
N<-1
B<-2/N
tru<-apply(rbind(x,y),2,naivelog)
for(rel in reliability){
  A<-B^2*(1-rel)*V/rel
  if(is.nan(A)){
    A<-0
    B<-1
  }
  MU<-A/B
  Skew<-2*A/B^3
  Va=A/B^2
  Kurt=3*A*(2+A)/B^4
  u<-matrix(rgamma(n*samp,A,B),nrow=n,ncol=samp)
  w<-(x+u)
  data<-rbind(w,y)
  nai<-apply(data,2,naivelog)
  mcs<-apply(data,2,mccs)
  cons<-tess(Va,Skew,Kurt)
  b<-cons[1]
  c<-cons[2]
  d<-cons[3]
  sdw<-cons[4]
  print(cons[4])
  print(err(b,c,d))
  new<-apply(data,2,newlog)
  suff<-apply(data,2,sufflog)
  DATA<-data
  newbias1<-append(newbias1,mean(new[1,]-b0))
  newmse1<-append(newmse1,mean((new[1,]-b0)^2))
  naibias1<-append(naibias1,mean(nai[1,]-b0))
  naimse1<-append(naimse1,mean((nai[1,]-b0)^2))
  newbias2<-append(newbias2,mean(new[2,]-b1))
  newmse2<-append(newmse2,mean((new[2,]-b1)^2))
  naibias2<-append(naibias2,mean(nai[2,]-b1))
  naimse2<-append(naimse2,mean((nai[2,]-b1)^2))
  mccsbias1<-append(mccsbias1,mean((mcs[1,]-b0)))
  mccsmse1<-append(mccsmse1,mean((mcs[1,]-b0)^2))
  mccsbias2<-append(mccsbias2,mean((mcs[2,]-b1)))
  mccsmse2<-append(mccsmse2,mean((mcs[2,]-b1)^2))
  suffbias1<-append(suffbias1,mean((suff[1,]-b0)))
  suffmse1<-append(suffmse1,mean((suff[1,]-b0)^2))
  suffbias2<-append(suffbias2,mean((suff[2,]-b1)))
  suffmse2<-append(suffmse2,mean((suff[2,]-b1)^2))
  newiqr1<-append(newiqr1,IQR(new[1,]))
  newiqr2<-append(newiqr2,IQR(new[2,]))
  naiiqr1<-append(naiiqr1,IQR(nai[1,]))
  naiiqr2<-append(naiiqr2,IQR(nai[2,]))
  mccsiqr1<-append(mccsiqr1,IQR(mcs[1,]))
  mccsiqr2<-append(mccsiqr2,IQR(mcs[2,]))
  suffiqr1<-append(suffiqr1,IQR(suff[1,]))
  suffiqr2<-append(suffiqr2,IQR(suff[2,]))
  print(naimse1)
  print(mccsmse1)
  print(suffmse1)
  print(newmse1)
  print(naimse2)
  print(mccsmse2)
  print(suffmse2)
  print(newmse2)
  print(rel)
}

truiqr1<-IQR(tru[1,])
truiqr2<-IQR(tru[2,])
par(mfrow=c(2,3))
par(mar=c(4,3,3,1))
plot(NULL,xlim=c(.8,1),ylim=c(min(c(naibias1,newbias1,mccsbias1)),max(c(naibias1,newbias1))),main=TeX(r'(Bias of $\hat{\beta}_{0}$)'),ylab="",xlab="Reliability",las=1)
lines(reliability,naibias1,col="blue",lty=2,lwd=1.5)
lines(reliability,newbias1,col="red",lty=1,lwd=1.5)
lines(reliability,mccsbias1,col="brown",lty=4,lwd=1.5)
lines(reliability,suffbias1,col="chartreuse4",lty=5,lwd=1.5)
lines(reliability,rep(mean(tru[1,]-b0),length(reliability)),col="orange",lty=3,lwd=1.5)
points(reliability,naibias1,pch=0)
points(reliability,newbias1,pch=1)
points(reliability,mccsbias1,pch=2)
points(reliability,rep(mean((tru[1,]-b0)),length(reliability)),pch=3)
points(reliability,suffbias1,pch=4)

plot(NULL,xlim=c(.8,1),ylim=c(min(newmse1,suffmse1,mccsmse1),4),main=TeX(r'(Mean Squared Error of $\hat{\beta}_{0}$)'),ylab="",xlab="Reliability",las=1)
lines(reliability,naimse1,col="blue",lty=2,lwd=1.5)
lines(reliability,newmse1,col="red",lty=1,lwd=1.5)
lines(reliability,mccsmse1,col="brown",lty=4,lwd=1.5)
lines(reliability,rep(mean((tru[1,]-b0)^2),length(reliability)),col="orange",lty=3,lwd=1.5)
lines(reliability,suffmse1,col="chartreuse4",lty=5,lwd=1.5)
points(reliability,naimse1,pch=0)
points(reliability,newmse1,pch=1)
points(reliability,mccsmse1,pch=2)
points(reliability,rep(mean((tru[1,]-b0)^2),length(reliability)),pch=3)
points(reliability,suffmse1,pch=4)

plot(NULL,xlim=c(.8,1),ylim=c(min(c(newiqr1,suffiqr1,mccsiqr1)),max(c(naiiqr1,newiqr1,mccsiqr1))),main=TeX(r'(IQR of $\hat{\beta}_{0}$)'),ylab="",xlab="Reliability",las=1)
lines(reliability,naiiqr1,col="blue",lty=2,lwd=1.5)
lines(reliability,newiqr1,col="red",lty=1,lwd=1.5)
lines(reliability,mccsiqr1,col="brown",lty=4,lwd=1.5)
lines(reliability,rep(truiqr1,length(reliability)),col="orange",lty=3,lwd=1.5)
lines(reliability,suffiqr1,col="chartreuse4",lty=5,lwd=1.5)
points(reliability,naiiqr1,pch=0)
points(reliability,newiqr1,pch=1)
points(reliability,mccsiqr1,pch=2)
points(reliability,rep(truiqr1,length(reliability)),pch=3)
points(reliability,suffiqr1,pch=4)

plot(NULL,xlim=c(.8,1),ylim=c(min(c(naibias2,newbias2,mccsbias2)),max(c(naibias2,newbias2))),main=TeX(r'(Bias of $\hat{\beta}_{1}$)'),ylab="",xlab="Reliability",las=1)
lines(reliability,naibias2,col="blue",lty=2,lwd=1.5)
lines(reliability,newbias2,col="red",lty=1,lwd=1.5)
lines(reliability,mccsbias2,col="brown",lty=4,lwd=1.5)
lines(reliability,suffbias2,col="chartreuse4",lty=5,lwd=1.5)
lines(reliability,rep(mean(tru[2,]-b1),length(reliability)),col="orange",lty=3,lwd=1.5)
points(reliability,naibias2,pch=0)
points(reliability,newbias2,pch=1)
points(reliability,mccsbias2,pch=2)
points(reliability,rep(mean((tru[2,]-b1)),length(reliability)),pch=3)
points(reliability,suffbias2,pch=4)

plot(NULL,xlim=c(.8,1),ylim=c(min(newmse2,suffmse2,mccsmse2),4),main=TeX(r'(Mean Squared Error of $\hat{\beta}_{1}$)'),ylab="",xlab="Reliability",las=1)
lines(reliability,naimse2,col="blue",lty=2,lwd=1.5)
lines(reliability,newmse2,col="red",lty=1,lwd=1.5)
lines(reliability,mccsmse2,col="brown",lty=4,lwd=1.5)
lines(reliability,rep(mean((tru[2,]-b1)^2),length(reliability)),col="orange",lty=3,lwd=1.5)
lines(reliability,suffmse2,col="chartreuse4",lty=5,lwd=1.5)
points(reliability,naimse2,pch=0)
points(reliability,newmse2,pch=1)
points(reliability,mccsmse2,pch=2)
points(reliability,rep(mean((tru[2,]-b1)^2),length(reliability)),pch=3)
points(reliability,suffmse2,pch=4)

plot(NULL,xlim=c(.8,1),ylim=c(min(c(newiqr2,suffiqr2,mccsiqr2)),max(c(naiiqr2,newiqr2,mccsiqr2))),main=TeX(r'(IQR of $\hat{\beta}_{1}$)'),ylab="",xlab="Reliability",las=1)
lines(reliability,naiiqr2,col="blue",lty=2,lwd=1.5)
lines(reliability,newiqr2,col="red",lty=1,lwd=1.5)
lines(reliability,mccsiqr2,col="brown",lty=4,lwd=1.5)
lines(reliability,rep(truiqr2,length(reliability)),col="orange",lty=3,lwd=1.5)
lines(reliability,suffiqr2,col="chartreuse4",lty=5,lwd=1.5)
points(reliability,naiiqr2,pch=0)
points(reliability,newiqr2,pch=1)
points(reliability,mccsiqr2,pch=2)
points(reliability,rep(truiqr2,length(reliability)),pch=3)
points(reliability,suffiqr2,pch=4)
