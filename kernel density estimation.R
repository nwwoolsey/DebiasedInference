library(pracma)
library(deconvolve)
library(expm)
library(LambertW)
library(latex2exp)

tess<-function(V,S,K){
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
tess2<-function(V,S,K){
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
    out<-optim(c(sqrt(V),.1,-.1),fn=fn,gr=gr,method="BFGS",V=V,S=S,K=K)$par
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
    out<-optim(c(sqrt(V),.1,-.1),fn=fn,gr=gr,method="BFGS",V=V,S=S,K=K)$par
  }
  return(c(out,1))
}
lapcf<-function(t){
  b<-sqrt(Va/2)
  return(1/(1+b^2*(t)^2))
}
kern<-function(t,h){
  if(t != 0){
    return(pi*t*csch(pi*t)/lapcf(t/h))
  }
  if(t==0){
    return(1)
  }
}
kdedecon<-function(w,x,h,f,Va){
  input<-outer(w,x,FUN="-")
  m<-5000
  gamma<-.2
  k<-seq(0,m-1,1)
  xk<-(k-m/2)*gamma
  FFFT<-function(f,m,gamma,h){
    alpha1<-sqrt(2*pi/m)
    alpha2<-alpha1
    sg<-(seq(0,m-1,1)-m/2)*alpha1
    tj<-(seq(0,m-1,1)-m/2)*alpha2
    out<-alpha1*(-1)^k*fft((-1)^k*sapply(sg,f,h=h))
    return(out)
  }
  L0<-Re(FFFT(f,m,gamma,h=h))
  index<-which(xk<max(input/h) & xk>min(input/h))
  l0<-L0[index]
  xk<-xk[index]
  closest<-function(x,L,xk){
    return(L[which(abs(xk-x)==min(abs(xk-x)))])
  }
  out<-1/(n*h)*matrix(vapply(input/h,closest,numeric(1),L=l0,xk=xk),ncol=dim(input)[2],nrow=dim(input)[1])
  out<-Re(apply(out,2,sum))
  out[which(out<0)]<-0
  out<-out/sum(out*.1)
  return(out)
}


kdedeconh<-function(data,x,hs,f){
  tol<-10^10
  for(h in hs){
    temp<-kdedecon(data,x,h,kern,Va) #deconvolve(data,xx=x,phiU=lapcf,bw=h,rescale=TRUE)$pdf##
    error<-mean((temp-true(x))^2)
    if(is.na(error)==TRUE){#if bandwidth is too small temp will be NaN
      next
    }
    if(error<tol){
      hopt<-h
      tol<-error
      tempopt<-temp
    }
  }
  out<-c()
  out$fit<-tempopt
  out$h<-hopt
  out$mise<-mean((tempopt-true(x))^2)
  return(c(out$mise,out$h))
}
kde<-function(data,x,h,b,c,d){
  K<-function(x,h,b,c,d){
    out<-Re(exp((b*1i+c+d*1i+x)/h)*(1+exp(2*(c+d*1i)/h)+4*exp((b*1i+c+d*1i+x)/h)+exp(2*(b*1i+x)/h)*(1+exp(2*(c+d*1i)/h)))/(2*(exp((c+d*1i)/h)+exp((b*1i+x)/h))^2*(1+exp((b*1i+c+d*1i+x)/h))^2))
    return(out)
  }
  input<-outer(data,x,FUN="-")
  out<-K(input,h,b,c,d)
  out[out<0]<-0
  out<-apply(out/n/h,2,sum)
  out<-out/sum(out*.1)
  return(out)
}
kdeh<-function(data,x,hs){
  tol<-10^10
  MU<-0
  xprime<-(x-MU)/scale
  dataprime<-(data-MU)/scale
  for(h in hs){
    temp<-kde(dataprime,xprime,h,b,c,d)
    error<-mean((temp-true(x))^2)
    if(is.na(error)==TRUE){#if bandwidth is too small temp will be NaN
      next
    }
    if(error<tol){
      hopt<-h
      tol<-error
      tempopt<-temp
    }
  }
  out<-c()
  out$fit<-tempopt
  out$h<-hopt
  out$mise<-tol
  return(c(out$mise,out$h))
}
kdenai<-function(data,x,h){
  #K<-function(x,h){
  #  out<-1/sqrt(2*pi)*exp(-x^2/2/h^2)
  #  return(out)
  #}
  K<-function(x,h){
    out<-1/(2+exp(x/h)+exp(-x/h))
    return(out)
  }
  inputs<-outer(data,x,FUN="-")
  out<-apply(K(inputs,h)/n/h,2,sum)
  out<-out/sum(out*.1)
  return(out)
}
kdenaih<-function(data,x,hs){
  tol<-10^10
  df<-approxfun(density(w-avg))
  tar<-df(x)
  tar[which(is.na(tar))]=0
  for(h in hs){
    temp<-kdenai(data,x,h)
    error<-mean((temp-tar)^2)
    if(is.na(error)==TRUE){#if bandwidth is too small temp will be NaN
      next
    }
    if(error<tol){
      hopt<-h
      tol<-error
      tempopt<-temp
    }
  }
  out<-c()
  out$fit<-tempopt
  out$h<-hopt
  out$mise<-mean((tempopt-true(x))^2)
  return(c(out$mise,out$h))
}

true<-function(x){
  return(.5*dnorm(x,0,1)+.5*dnorm(x,5,1))
}


kdetru<-function(data,x,h){
  K<-function(x){
    out<-1/(2+exp(x/h)+exp(-x/h))
    #out<-1/sqrt(2*pi)*exp(-x^2/2/h^2)
    return(out)
  }
  inputs<-outer(data,x,FUN="-")
  out<-apply(K(inputs)/n/h,2,sum)
  out<-out/sum(out*.1)
  return(out)
}
kdetruh<-function(data,x,hs){
  tol<-10^10
  for(h in hs){
    temp<-kdetru(data,x,h)
    error<-mean((temp-true(x))^2)
    if(is.na(error)==TRUE){#if bandwidth is too small temp will be NaN
      next
    }
    if(error<tol){
      hopt<-h
      tol<-error
      tempopt<-temp
    }
  }
  out<-c()
  out$fit<-tempopt
  out$h<-hopt
  out$mise<-mean((tempopt-true(x))^2)
  return(c(out$mise,out$h))
}

par(mfrow=c(2,3))
par(mar=c(4,3,3,1))
set.seed(1)
samp<-500
n<-50
p<-matrix(rbinom(n*samp,1,.5),nrow=n,ncol=samp)
x<-matrix(p*rnorm(n*samp,0,1)+(1-p)*rnorm(n*samp,5,1),nrow=n,ncol=samp)
inputs<-seq(-2,7,9/99)
reliability<-c(.75,.85,.95,1)
newmse<-c()
naimse<-c()
trumse<-c()
decmse<-c()
new2mse<-c()
#compmse<-c()
for(rel in reliability){
  Vau<-var(matrix(x,ncol=1))[1]*(1-rel)/rel
  #s<-log(-1+Vau/(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)+(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)/Vau)#log(1/2*exp(-2*mu)*(exp(2*mu)+e xp(mu)*sqrt(exp(2*mu)+4*(1/rel-1)*var(matrix(x,ncol=1))[1])))
  #mu<--log(4)+log((sqrt(-1+Vau/(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)+(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)/Vau)*(Vau^4+2*sqrt(Vau^4*(4+Vau))*(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)+Vau^3*(4+(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3))+Vau^2*(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)*(4+(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3))+2*Vau*(sqrt(Vau^4*(4+Vau))+2*(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(2/3))))/((4+Vau)*(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(2/3)))
  temp<-2^(1/3)*Vau/((Vau^2*(N^2+2*Vau)+sqrt(N^2*Vau^4*(N^2+4*Vau)))^(1/3))
  sigma<-sqrt(2)*sqrt(log(sqrt(-1+temp+1/temp)))
  s<-sigma^2
  mu<-1/2*(-sigma^2+log(Vau/(exp(sigma^2)-1)))
  if(is.nan(s)){
    s<-0
  }
  if(is.nan(mu)){
    mu<-0
  }
  u<-matrix(rlnorm(n*samp,mu,sqrt(s)),nrow=n,ncol=samp)
  avg<-exp(mu+s/2)
  Skew<-exp(3*mu+3*s/2)*(exp(s)-1)^2*(2+exp(s))
  Va=exp(2*mu+s)*(exp(s)-1)
  Kurt=exp(4*mu+2*s)*(exp(s)-1)^2*(-3+exp(2*s)*(3+exp(s)*(2+exp(s))))
  print(c(Va,Skew,Kurt))
  o<-tess(Va,Skew,Kurt)
  scale<-o[4]
  b<-o[1]
  c<-o[2]
  d<-o[3]
  w<-x+u
  dec<-future.apply::future_apply(w-avg,2,kdedeconh,x=inputs,h=seq(.05,.45,.01),f=kern)
  print(mean(dec[1,]))
  new<-future.apply::future_apply(w-avg,2,kdeh,x=inputs,hs=seq(.02,.8,.01))
  o<-tess2(Va,Skew,Kurt)
  b<-o[1]
  c<-o[2]
  d<-o[3]
  scale<-1
  new2<-future.apply::future_apply(w-avg,2,kdeh,x=inputs,hs=seq(.02,.8,.01))
  print(mean(new[1,]))
  #com<-future.apply::future_apply(w-avg,2,kdehcomp,x=inputs,hs=seq(.02,.6,.01))
  print(mean(com[1,]))
  nai<-future.apply::future_apply(w-avg,2,kdenaih,x=inputs,hs=seq(.1,7,.1))
  print(mean(nai[1,]))
  newmse<-append(newmse,mean(new[1,]))
  naimse<-append(naimse,mean(nai[1,]))
  decmse<-append(decmse,mean(dec[1,]))
  new2mse<-append(new2mse,mean(new2[1,]))
  #compmse<-append(compmse,mean(com[1,]))
  mse<-cbind(newmse,naimse,decmse,new2mse)
  write.csv(x=mse,file="mse.csv")
  #print(rel)
}
tru<-apply(x,2,kdetruh,x=inputs,hs=seq(.1,5,.1))
trumse<-rep(mean(tru[1,]),length(reliability))

plot(NULL,xlim=c(.75,1),ylim=c(min(c(trumse)),max(mse[,2])),xlab="Reliability",ylab="",main=paste("MISE n=",n,sep=" "),las=1)
lines(reliability,mse[,1],col="red",lty=1,lwd=2)
lines(reliability,mse[,2],col="blue",lty=2,lwd=2)
lines(reliability,mse[,3],col="chartreuse4",lty=4,lwd=2)
lines(reliability,trumse,col="orange",lty=3,lwd=2)



n<-75
p<-matrix(rbinom(n*samp,1,.5),nrow=n,ncol=samp)
x<-matrix(p*rnorm(n*samp,0,1)+(1-p)*rnorm(n*samp,5,1),nrow=n,ncol=samp)
inputs<-seq(-2,8,.1)
reliability<-c(.75,.85,.95,1)
newmse2<-c()
naimse2<-c()
trumse2<-c()
decmse2<-c()
new2mse2<-c()
#compmse<-c()
for(rel in reliability){
  Vau<-var(matrix(x,ncol=1))[1]*(1-rel)/rel
  #s<-log(-1+Vau/(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)+(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)/Vau)#log(1/2*exp(-2*mu)*(exp(2*mu)+e xp(mu)*sqrt(exp(2*mu)+4*(1/rel-1)*var(matrix(x,ncol=1))[1])))
  #mu<--log(4)+log((sqrt(-1+Vau/(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)+(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)/Vau)*(Vau^4+2*sqrt(Vau^4*(4+Vau))*(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)+Vau^3*(4+(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3))+Vau^2*(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)*(4+(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3))+2*Vau*(sqrt(Vau^4*(4+Vau))+2*(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(2/3))))/((4+Vau)*(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(2/3)))
  temp<-2^(1/3)*Vau/((Vau^2*(N^2+2*Vau)+sqrt(N^2*Vau^4*(N^2+4*Vau)))^(1/3))
  sigma<-sqrt(2)*sqrt(log(sqrt(-1+temp+1/temp)))
  s<-sigma^2
  mu<-1/2*(-sigma^2+log(Vau/(exp(sigma^2)-1)))
  if(is.nan(s)){
    s<-0
  }
  if(is.nan(mu)){
    mu<-0
  }
  u<-matrix(rlnorm(n*samp,mu,sqrt(s)),nrow=n,ncol=samp)
  avg<-exp(mu+s/2)
  Skew<-exp(3*mu+3*s/2)*(exp(s)-1)^2*(2+exp(s))
  Va=exp(2*mu+s)*(exp(s)-1)
  Kurt=exp(4*mu+2*s)*(exp(s)-1)^2*(-3+exp(2*s)*(3+exp(s)*(2+exp(s))))
  print(c(Va,Skew,Kurt))
  o<-tess(Va,Skew,Kurt)
  scale<-o[4]
  b<-o[1]
  c<-o[2]
  d<-o[3]
  w<-x+u
  dec<-future.apply::future_apply(w-avg,2,kdedeconh,x=inputs,h=seq(.05,.45,.01),f=kern)
  print(mean(dec[1,]))
  new<-future.apply::future_apply(w-avg,2,kdeh,x=inputs,hs=seq(.02,.8,.01))
  o<-tess2(Va,Skew,Kurt)
  b<-o[1]
  c<-o[2]
  d<-o[3]
  scale<-1
  new2<-future.apply::future_apply(w-avg,2,kdeh,x=inputs,hs=seq(.02,.8,.01))
  print(mean(new[1,]))
  #com<-future.apply::future_apply(w-avg,2,kdehcomp,x=inputs,hs=seq(.02,.6,.01))
  print(mean(com[1,]))
  nai<-future.apply::future_apply(w-avg,2,kdenaih,x=inputs,hs=seq(.1,7,.1))
  print(mean(nai[1,]))
  newmse2<-append(newmse,mean(new[1,]))
  naimse2<-append(naimse,mean(nai[1,]))
  decmse2<-append(decmse,mean(dec[1,]))
  new2mse2<-append(new2mse,mean(new2[1,]))
  #compmse<-append(compmse,mean(com[1,]))
  mse2<-cbind(newmse,naimse,decmse,new2mse)
  write.csv(x=mse,file="mse2.csv")
  #print(rel)
}
tru<-apply(x,2,kdetruh,x=inputs,hs=seq(.1,5,.1))
trumse<-rep(mean(tru[1,]),length(reliability))

plot(NULL,xlim=c(.75,1),ylim=c(min(c(trumse)),max(mse[,2])),xlab="Reliability",ylab="",main=paste("MISE n=",n,sep=" "),las=1)
lines(reliability,mse[,1],col="red",lty=1,lwd=2)
lines(reliability,mse[,2],col="blue",lty=2,lwd=2)
lines(reliability,mse[,3],col="chartreuse4",lty=4,lwd=2)
lines(reliability,trumse,col="orange",lty=3,lwd=2)


#############################################################################################
  
  kdedecon<-function(w,x,h,f,Va){
    input<-outer(w,x,FUN="-")
    m<-5000
    gamma<-.2
    k<-seq(0,m-1,1)
    xk<-(k-m/2)*gamma
    FFFT<-function(f,m,gamma,h){
      alpha1<-sqrt(2*pi/m)
      alpha2<-alpha1
      sg<-(seq(0,m-1,1)-m/2)*alpha1
      tj<-(seq(0,m-1,1)-m/2)*alpha2
      out<-alpha1*(-1)^k*fft((-1)^k*sapply(sg,f,h=h))
      return(out)
    }
    L0<-Re(FFFT(f,m,gamma,h=h))
    index<-which(xk<max(input/h) & xk>min(input/h))
    l0<-L0[index]
    xk<-xk[index]
    closest<-function(x,L,xk){
      return(L[which(abs(xk-x)==min(abs(xk-x)))])
    }
    out<-1/(n*h)*matrix(vapply(input/h,closest,numeric(1),L=l0,xk=xk),ncol=dim(input)[2],nrow=dim(input)[1])
    out<-Re(apply(out,2,sum))
    out[which(out<0)]<-0
    out<-out/sum(out*.1)
    return(out)
  }


kdedeconh<-function(data,x,hs,f){
  tol<-10^10
  for(h in hs){
    temp<-kdedecon(data,x,h,kern,Va) #deconvolve(data,xx=x,phiU=lapcf,bw=h,rescale=TRUE)$pdf##
    error<-mean((temp-true(x))^2)
    if(is.na(error)==TRUE){#if bandwidth is too small temp will be NaN
      next
    }
    if(error<tol){
      hopt<-h
      tol<-error
      tempopt<-temp
    }
  }
  out<-c()
  out$fit<-tempopt
  out$h<-hopt
  out$mise<-mean((tempopt-true(x))^2)
  return(out$fit)
}
kde<-function(data,x,h,b,c,d){
  #K<-function(x,h,b,c,d){
  #  out<-1/sqrt(2*pi)*1/2*(exp(1/2*(-(x-c)^2/h^2+(b-d)^2/h^2))*cos((x-c)*(b-d)/h^2)+exp(1/2*(-(x+c)^2/h^2+(b+d)^2/h^2))*cos((x+c)*(b+d)/h^2))
  #  return(out)
  #}
  K<-function(x,h,b,c,d){
    #out<-Re(exp((b*1i+c+d*1i+x)/h)*(1+exp(2*(c+d*1i)/h)+4*exp((b*1i+c+d*1i+x)/h)+exp(2*(b*1i+x)/h)*(1+exp(2*(c+d*1i)/h)))/(2*(exp((c+d*1i)/h)+exp((b*1i+x)/h))^2*(1+exp((b*1i+c+d*1i+x)/h))^2))
    out<-Re(1/2*(exp(-1/2*(x+1i*(b+1i*c-d))^2)+exp(-1/2*(x+b*1i+c+d*1i)^2)))
    return(out)
  }
  input<-outer(data,x,FUN="-")
  out<-K(input,h,b,c,d)
  out[out<0]<-0
  out<-apply(out/n/h,2,sum)
  out<-out/sum(out*.1)
  return(out)
}
kdeh<-function(data,x,hs){
  tol<-10^10
  MU<-0
  xprime<-(x-MU)/scale
  dataprime<-(data-MU)/scale
  for(h in hs){
    temp<-kde(dataprime,xprime,h,b,c,d)
    error<-mean((temp-true(x))^2)
    if(is.na(error)==TRUE){#if bandwidth is too small temp will be NaN
      next
    }
    if(error<tol){
      hopt<-h
      tol<-error
      tempopt<-temp
    }
  }
  out<-c()
  out$fit<-tempopt
  out$h<-hopt
  out$mise<-tol
  return(out$fit)
}
kdenai<-function(data,x,h){
  #K<-function(x,h){
  #  out<-1/sqrt(2*pi)*exp(-x^2/2/h^2)
  #  return(out)
  #}
  K<-function(x,h){
    out<-1/(2+exp(x/h)+exp(-x/h))
    return(out)
  }
  inputs<-outer(data,x,FUN="-")
  out<-apply(K(inputs,h)/n/h,2,sum)
  out<-out/sum(out*.1)
  return(out)
}
kdenaih<-function(data,x,hs){
  tol<-10^10
  df<-approxfun(density(w-avg))
  tar<-df(x)
  tar[which(is.na(tar))]=0
  for(h in hs){
    temp<-kdenai(data,x,h)
    error<-mean((temp-tar)^2)
    if(is.na(error)==TRUE){#if bandwidth is too small temp will be NaN
      next
    }
    if(error<tol){
      hopt<-h
      tol<-error
      tempopt<-temp
    }
  }
  out<-c()
  out$fit<-tempopt
  out$h<-hopt
  out$mise<-mean((tempopt-true(x))^2)
  return(out$fit)
}

true<-function(x){
  return(.5*dnorm(x,0,1)+.5*dnorm(x,5,1))
}


kdetru<-function(data,x,h){
  K<-function(x){
    out<-1/(2+exp(x/h)+exp(-x/h))
    #out<-1/sqrt(2*pi)*exp(-x^2/2/h^2)
    return(out)
  }
  inputs<-outer(data,x,FUN="-")
  out<-apply(K(inputs)/n/h,2,sum)
  out<-out/sum(out*.1)
  return(out)
}
kdetruh<-function(data,x,hs){
  tol<-10^10
  for(h in hs){
    temp<-kdetru(data,x,h)
    error<-mean((temp-true(x))^2)
    if(is.na(error)==TRUE){#if bandwidth is too small temp will be NaN
      next
    }
    if(error<tol){
      hopt<-h
      tol<-error
      tempopt<-temp
    }
  }
  out<-c()
  out$fit<-tempopt
  out$h<-hopt
  out$mise<-mean((tempopt-true(x))^2)
  return(out$fit)
}
n<-50
p<-matrix(rbinom(n*samp,1,.5),nrow=n,ncol=samp)
x<-matrix(p*rnorm(n*samp,0,1)+(1-p)*rnorm(n*samp,5,1),nrow=n,ncol=samp)
inputs<-seq(-2,8,.1)
reliability<-c(.75,.85,.95,1)
newmse<-c()
naimse<-c()
trumse<-c()
decmse<-c()
new2mse<-c()
#compmse<-c()
rel=.9
Vau<-var(matrix(x,ncol=1))[1]*(1-rel)/rel
#s<-log(-1+Vau/(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)+(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)/Vau)#log(1/2*exp(-2*mu)*(exp(2*mu)+e xp(mu)*sqrt(exp(2*mu)+4*(1/rel-1)*var(matrix(x,ncol=1))[1])))
#mu<--log(4)+log((sqrt(-1+Vau/(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)+(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)/Vau)*(Vau^4+2*sqrt(Vau^4*(4+Vau))*(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)+Vau^3*(4+(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3))+Vau^2*(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3)*(4+(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(1/3))+2*Vau*(sqrt(Vau^4*(4+Vau))+2*(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(2/3))))/((4+Vau)*(4*sqrt(Vau^4*(4+Vau))+Vau^2*(8+Vau))^(2/3)))
temp<-2^(1/3)*Vau/((Vau^2*(N^2+2*Vau)+sqrt(N^2*Vau^4*(N^2+4*Vau)))^(1/3))
sigma<-sqrt(2)*sqrt(log(sqrt(-1+temp+1/temp)))
s<-sigma^2
mu<-1/2*(-sigma^2+log(Vau/(exp(sigma^2)-1)))
if(is.nan(s)){
  s<-0
}
if(is.nan(mu)){
  mu<-0
}
u<-matrix(rlnorm(n*samp,mu,sqrt(s)),nrow=n,ncol=samp)
avg<-exp(mu+s/2)
Skew<-exp(3*mu+3*s/2)*(exp(s)-1)^2*(2+exp(s))
Va=exp(2*mu+s)*(exp(s)-1)
Kurt=exp(4*mu+2*s)*(exp(s)-1)^2*(-3+exp(2*s)*(3+exp(s)*(2+exp(s))))
print(c(Va,Skew,Kurt))
o<-tess(Va,Skew,Kurt)
scale<-o[4]
b<-o[1]
c<-o[2]
d<-o[3]
w<-x+u
#new<-future.apply::future_apply(w-avg,2,kdeh,x=inputs,hs=seq(.02,.8,.01))
dec<-future.apply::future_apply(w-avg,2,kdedeconh,x=inputs,h=seq(.05,.45,.01),f=kern)
print(mean(dec[1,]))
o<-tess2(Va,Skew,Kurt)
b<-o[1]
c<-o[2]
d<-o[3]
scale<-1
new2<-future.apply::future_apply(w-avg,2,kdeh,x=inputs,hs=seq(.02,.8,.01))
nai<-future.apply::future_apply(w-avg,2,kdenaih,x=inputs,hs=seq(.1,7,.1))
tru<-future.apply::future_apply(x,2,kdetruh,x=inputs,hs=seq(.01,2,.1))
boxplot(t(dec),pch=".",ylim=c(0,.3),xaxt="n",yaxt="n",main="Deconvoluting Kernel")
lines(seq(1,91,.01),true(seq(1,91,.01)*.1-2.12),lwd=2.0)
boxplot(t(new2),pch=".",ylim=c(0,.3),xaxt="n",yaxt="n",main="Tessarine Correction")
lines(seq(1,91,.01),true(seq(1,91,.01)*.1-2.12),lwd=2.0)
boxplot(t(nai),pch=".",ylim=c(0,.3),xaxt="n",yaxt="n",main="Naive Estimator")
lines(seq(1,91,.01),true(seq(1,91,.01)*.1-2.12),lwd=2.0)
boxplot(t(tru),pch=".",ylim=c(0,.3),xaxt="n",yaxt="n",main="True Data Estimator")
lines(seq(1,91,.01),true(seq(1,91,.01)*.1-2.12),lwd=2.0)

