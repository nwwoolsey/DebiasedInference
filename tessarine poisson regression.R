###File to generate Poisson regression plot

library(latex2exp)
tess<-function(V,S,K){#calculates tessarine coefficients
  eq1<-c(0,0,0,0,S^2,0,45*V^2-9*K,0,-36*V,0,36)
  if(V==0){
    return(c(0,0,0,0))
  }
  sols<-polyroot(eq1)
  ds<- Re(sols[which(abs(Im(sols))<10^(-4))])
  if(length(ds)==0){
    break
  }
  min<-1000000000000000
  for(d in ds){
    eq2<-c(-S^2,0,36*d^2*(V-d^2),0,36*d^2)
    sols2<-polyroot(eq2)
    cs<- Re(sols2[which(abs(Im(sols2))<10^(-4))])
    if(length(cs)==0){
      break
    }
    for(c in cs){
      b<-sqrt(V+c^2-d^2)
      temp<-sqrt(b^2+c^2+d^2)
      if(temp< min){
        out<-c(b,c,d)
        min<-temp
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
    f1<-function(b,c,d){
      return(V-b^2+c^2-d^2)
    }
    f2<-function(b,c,d){
      return(S-6*b*c*d)
    }
    f3<-function(b,c,d){
      return(K+b^4+c^4+d^4-6*V*(b^2-c^2+d^2)-6*(b^2*c^2-b^2*d^2+c^2*d^2))
    }
    b<-out[1]
    c<-out[2]
    d<-out[3]
    err<-abs(f1(b,c,d))/V+abs(f2(b,c,d)/3/S)+abs(f3(b,c,d)/6)/K
  return(c(out,err))
}

naivepois<-function(data){#naive Poisson Regression
  x<-data[1:n]-MU
  y<-data[(n+1):(2*n)]
  score<-function(theta){
   b0<-theta[1]
    b1<-theta[2]
    score1<--sum(y-exp(b0+b1*x))
    score2<--sum(x*y-exp(b0+b1*x)*x)
    return(score1^2+score2^2)
  }
  out<-optim(c(0,0),fn=score,method="BFGS")$par
  return(out)
}
naivepois<-function(data){#naive estimator
  x<-data[1:n]-mu
  y<-data[(n+1):(2*n)]
  xbar<-mean(x)
  sdx<-sd(x)
  xstar<-(x-xbar)/sdx
  likelihood<-function(theta){
    b0<-theta[1]
    b1<-theta[2]
    return(-sum(y*(b0+b1*xstar)-exp(b0+b1*xstar)))
  }
  score<-function(theta){
    b0<-theta[1]
    b1<-theta[2]
    score1<--sum(y-exp(b0+b1*xstar))
    score2<--sum(xstar*y-exp(b0+b1*xstar)*xstar)
    return(c(score1,score2))
  }
  out<-optim(c(0,0),fn=likelihood,gr=score,method="BFGS")$par
  b0<-out[1]-out[2]*xbar/sdx
  b1<-out[2]/sdx
  return(c(b0,b1))
}
pois<-function(data){#Novel Poisson Regression
  w<-data[1:n]-MU
  y<-data[(n+1):(2*n)]
  muw<-0#-min(w)#mean(w)
  sdw<-1#max(w)-min(w)
  Vastar<-Va/sdw^2
  Skstar<-Skew/sdw^3
  Kustar<-Kurt/sdw^4
  cons<-tess(Vastar,Skstar,Kustar)
  wstar<-(w-muw)/sdw
  b<-cons[1]
  c<-cons[2]
  d<-cons[3]
  score<-function(theta){
    b0<-theta[1]
    b1<-theta[2]
    score1<-mean(y-Re(exp(b0+b1*(wstar+b*1i))*cosh(b1*(c+d*1i))))^2
    score2<-mean(y*wstar-Re(exp(b0+b1*(wstar+b*1i))*((wstar+b*1i)*cosh(b1*(c+d*1i))+(c+d*1i)*sinh(b1*(c+d*1i)))))^2
   return(score1+score2)
  }
  nai<-naivepois(data)
  out<-optim(nai,fn=score,method="BFGS")$par
  b0<-out[1]
  b1<-out[2]
  return(c(b0,b1,b,c,d))
}


cs<-function(data){#corrected score
  w<-data[1:n]-MU
  y<-data[(n+1):(2*n)]
  muw<-0#mean(w)
  sdw<-1#sqrt(var(w))
  wstar<-(w-muw)/sdw
  Vastar<-Va/sdw^2
  score<-function(theta){
    b0<-theta[1]
    b1<-theta[2]
    score0<-mean(y-exp(b0-b1^2/2*Vastar+b1*wstar))
    score1<-mean(wstar*y-exp(b0-Vastar*b1^2/2+b1*wstar)*(wstar-Vastar*b1))
    return(score0^2+score1^2)
  }
  nai<-naivepois(data)
  out<-optim(nai,fn=score,method="BFGS")$par
  out[1]<-out[1]-out[2]*muw/sdw
  out[2]<-out[2]/sdw
  return(out)
}
n<-500
samp<-1000
set.seed(1)
reliability<-c(seq(.8,1,.1))
newbias0<-c()
newmse0<-c()
naivebias0<-c()
naivemse0<-c()
truebias0<-c()
truemse0<-c()
newbias1<-c()
newmse1<-c()
naivebias1<-c()
naivemse1<-c()
truebias1<-c()
truemse1<-c()
stefbias1<-c()
stefbias0<-c()
stefmse1<-c()
stefmse0<-c()
newiqr0<-c()
newiqr1<-c()
trueiqr0<-c()
trueiqr1<-c()
naiveiqr0<-c()
naiveiqr1<-c()
stefiqr0<-c()
stefiqr1<-c()


b0<-1
b1<--1
x<-matrix(rnorm(n*samp,0,1),nrow=n,ncol=samp)#True data 
lambda<-exp(b0+b1*x)#average Poisson responses
y<-apply(lambda,2,rpois,n=n)#Poisson responses
N<-1
B<-2/N
#B<-2
mom<-c()
con<-c()
shrink<-c()
V<-var(matrix(x,ncol=1))[1]#variance of true values
for(rel in reliability){
  A<-B^2*(1-rel)*V/rel
  u<-matrix(rgamma(n*samp,A,B),nrow=n,ncol=samp)#gamma errors
  w<-x+u
  MU<-A/B
  Va<-A/B^2
  Skew<-2*A/B^3
  Kurt<-3*A*(2+A)/B^4
  l<-tess(Va,Skew,Kurt)
  b<-l[1]
  c<-l[2]
  d<-l[3]
  ERR<-l[4]
  data<-rbind(w,y)#matrix of all data (covariates and responses)
  new<-apply(data,2,pois)#New method
  stef<-apply(data,2,cs)#corrected score
  nai<-apply(data,2,naivepois)#naive
  true<-apply(rbind(x+MU,y),2,naivepois)#True data estimate
  truebias1<-append(truebias1,(mean(true[2,]-b1)))
  truemse1<-append(truemse1,mean((true[2,]-b1)^2))
  newbias1<-append(newbias1,(mean(new[2,]-b1)))
  newmse1<-append(newmse1,mean((new[2,]-b1)^2))
  naivebias1<-append(naivebias1,(mean(nai[2,]-b1)))
  naivemse1<-append(naivemse1,mean((nai[2,]-b1)^2))
  truebias0<-append(truebias0,(mean((true[1,]-b0))))
  truemse0<-append(truemse0,mean((true[1,]-b0)^2))
  newbias0<-append(newbias0,(mean(new[1,]-b0)))
  newmse0<-append(newmse0,mean((new[1,]-b0)^2))
  naivebias0<-append(naivebias0,(mean(nai[1,]-b0)))
  naivemse0<-append(naivemse0,mean((nai[1,]-b0)^2))
  stefbias0<-append(stefbias0,(mean(stef[1,]-b0)))
  stefmse0<-append(stefmse0,mean((stef[1,]-b0)^2))
  stefbias1<-append(stefbias1,(mean(stef[2,]-b1)))
  stefmse1<-append(stefmse1,mean((stef[2,]-b1)^2))
  newiqr0<-append(newiqr0,IQR(new[1,]))
  newiqr1<-append(newiqr1,IQR(new[2,]))
  trueiqr0<-append(trueiqr0,IQR(true[1,]))
  trueiqr1<-append(trueiqr1,IQR(true[2,]))
  naiveiqr0<-append(naiveiqr0,IQR(nai[1,]))
  naiveiqr1<-append(naiveiqr1,IQR(nai[2,]))
  stefiqr0<-append(stefiqr0,IQR(stef[1,]))
  stefiqr1<-append(stefiqr1,IQR(stef[2,]))
  print(rel)
  print(newmse1)
  print(naivemse1)
  print(stefmse1)
}



par(mfrow=c(2,3))
par(mar=c(4,4,3,1))
plot(NULL,xlim=c(.8,1),ylim=c(-.1,max(naivebias0)),xlab="Reliability",ylab="",main=TeX(r'(Bias of $\hat{\beta}_{0}$)'),las=1)
lines(reliability,naivebias0,lty=2,col="blue",lwd=1.5)
lines(reliability,truebias0,lty=3,col="orange",lwd=1.5)
lines(reliability,stefbias0,lty=4,col="chartreuse4",lwd=1.5)
lines(reliability,newbias0,lty=1,col="red",lwd=1.5)



plot(NULL,xlim=c(.8,1),ylim=c(min(truemse0),.01),xlab="Reliability",ylab="",main=TeX(r'(Mean Squared Error of $\hat{\beta}_{0}$)'),las=1)
lines(reliability,naivemse0,lty=2,col="blue",lwd=1.5)
lines(reliability,newmse0,lty=1,col="red",lwd=1.5)
lines(reliability,truemse0,lty=3,col="orange",lwd=1.5)
lines(reliability,stefmse0,lty=4,col="chartreuse4",lwd=1.5)

plot(NULL,xlim=c(.8,1),ylim=c(.04,max(c(newiqr0,stefiqr0,naiveiqr0))),xlab="Reliability",ylab="", main=TeX(r'(IQR of $\hat{\beta}_{0}$)'),las=1)
lines(reliability,naiveiqr0,lty=2,col="blue",lwd=1.5)
lines(reliability,trueiqr0,lty=3,col="orange",lwd=1.5)
lines(reliability,stefiqr0,lty=4,col="chartreuse4",lwd=1.5)
lines(reliability,newiqr0,lty=1,col="red",lwd=1.5)


plot(NULL,xlim=c(.8,1),ylim=c(-.15,max(naivebias1)),xlab="Reliability",ylab="",main=TeX(r'(Bias of $\hat{\beta}_{1}$)'),las=1)
lines(reliability,naivebias1,lty=2,col="blue",lwd=1.5)
lines(reliability,truebias1,lty=3,col="orange",lwd=1.5)
lines(reliability,stefbias1,lty=4,col="chartreuse4",lwd=1.5)
lines(reliability,newbias1,lty=1,col="red",lwd=1.5)

plot(NULL,xlim=c(.8,1),ylim=c(min(truemse1),max(naivemse1)),xlab="Reliability",ylab="",main=TeX(r'(Mean Squared Error of $\hat{\beta}_{1}$)'),las=1)
lines(reliability,naivemse1,lty=2,col="blue",lwd=1.5)
lines(reliability,newmse1,lty=1,col="red",lwd=1.5)
lines(reliability,truemse1,lty=3,col="orange",lwd=1.5)
lines(reliability,stefmse1,lty=4,col="chartreuse4",lwd=1.5)


plot(NULL,xlim=c(.8,1),ylim=c(.024,max(c(newiqr1,stefiqr1,naiveiqr1))),xlab="Reliability",ylab="", main=TeX(r'(IQR of $\hat{\beta}_{1}$)'),las=1)
lines(reliability,naiveiqr1,lty=2,col="blue",lwd=1.5)
lines(reliability,trueiqr1,lty=3,col="orange",lwd=1.5)
lines(reliability,stefiqr1,lty=4,col="chartreuse4",lwd=1.5)
lines(reliability,newiqr1,lty=1,col="red",lwd=1.5)
