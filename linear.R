rm(list=ls())
ancdf=function(x,mu,sigma,tao){
  if (x <= mu) {
    z =(2/sqrt(sigma*pi))*((sqrt(1/(1-tao))+sqrt(1/tao))^(-1))*exp(-((1-tao)*(x-mu)^2)/sigma)
  }else {
    z =(2/sqrt(sigma*pi))*((sqrt(1/(1-tao))+sqrt(1/tao))^(-1))*exp(-((tao)*(x-mu)^2)/sigma)
  }
  return(z)
}
negLogLikensumli=function(y,x,beta,b.l,sigma,C.im,w.k){
  
  J <- dim(as.array(unique(y)))[1]
  n <- dim(x)[1]
  C.im=as.matrix(C.im)
  if(is.null(dim(C.im))){
    L=1
  }else{ L= dim(C.im)[2]}
  lnpdf1 <- array(0, dim = c(n, L))
  mu <- x %*% beta
  tau=rep(0,L)
  for(i.ta in 1:L){
    tau[i.ta]=i.ta/(1+L)
  }
  for (i in 1:n) {  
    for(j in 1:L){
      
      meanp <- mu[i]+b.l[j]
      an=ancdf(y[i], meanp, sigma, tau[j])
      if(an==0){
        lnpdf1[i,j]=0
      }
      if(an!=0){
        lnpdf1[i,j] <- C.im[i,j]*log(w.k[j]*ancdf(y[i], meanp, sigma, tau[j])) ## (x,mu,sigma,tao) 
      }
    }
  }
  lnpdf=rowSums(lnpdf1)
  nlogl <- -lnpdf
  negsumlogl <- -sum(lnpdf)
  respon <- list(nlogl = nlogl, negsumlogl = negsumlogl)
  return(respon)
}

#####DIC
dicli=function (y, x, betadraws, postMeanbeta, b.n,
                burn, mcmc,sigma,C.i,w.n) 
{
  cols <- colnames(x)
  names(x) <- NULL
  names(y) <- NULL
  x <- as.matrix(x)
  y <- as.matrix(y)
  nsim <-   mcmc-burn
  L=dim(b.n)[1]
  postb.l=apply(b.n,1,mean)
  betadraws <- as.matrix(betadraws)
  postsigma=mean(sigma)
  sumC.im=matrix(0,n,L)
  for(ij in 1:nsim){
    ij1=(1+(ij-1)*L):(ij*L)
    sumC.im=sumC.im+C.i[,ij1]
  }
  postC.im=sumC.im/nsim
  postw.n=apply(w.n,1,mean)
  dev <- array(0, dim = c(1))
  DIC <- array(0, dim = c(1))
  pd <- array(0, dim = c(1))
  ans <-  negLogLikensumli(y, x, postMeanbeta,postb.l , 
                           postsigma,postC.im,postw.n)
  dev <- 2 * ans$negsumlogl 
  Deviance <- array(0, dim = c(nsim, 1))
  for (i in 1:nsim) {
    
    ii=(1+(i-1)*L):(i*L)
    temp <-  negLogLikensumli(y, x, betadraws[,i],b.n[,i], sigma[i],C.i[,ii],w.n[,i])
    Deviance[i, 1] <- 2 * temp$negsumlogl
  }
  avgdDeviance <- mean(Deviance)
  DIC <- 2 * avgdDeviance - dev
  pd <- avgdDeviance - dev
  result <- list(DIC = DIC, pd = pd, dev = dev)
  return(result)
}


library(MASS) 
library(quantreg)
library(LaplacesDemon)
library(mvtnorm)
library(rmutil)
library(Brq)
library(spdep)
tt.g=rep(0,15)
ct.g=rep(0,15)
l.v=0
z.v=0
c.lv=0
c.zv=0
set.seed(202307)
p.x=1
RMSE.c=rep(0,p.x)
AE.c=rep(0,p.x)
AD.c=rep(0,p.x)
ME.c=rep(0,p.x)
beta.mm=matrix(0,15,p.x)
system.time({for(i.x in 1:p.x){

  data(boston, package="spData")
  y=as.matrix(log(boston.c$CMEDV))
  x.1=boston.c[,-c(1:3)]
  x.2=(x.1[,-c(3:4)])
  names(x.2)=NULL
  x=apply(x.2,2,as.numeric)
  x=scale(x,center=TRUE,scale=TRUE)
  y=scale(y,center=TRUE,scale=TRUE)
  p=ncol(x)
  n=length(y)
  ###########tau
  L=1
  tau=rep(0,L)
  for(i.ta in 1:L){
    tau[i.ta]=i.ta/(1+L)
  }
  ##########
  beta.c=rep(1,p)
  b.c=rep(1,L)
  sigma.c=1
  gama.c=rep(1,p)
  delta.c=rep(1,p)
  C.im=matrix(0,n,L)
  p.m1=rep(1/L,L)
  for(i.cc in 1:n){
    C.im[i.cc,]=rmultinom(1,1,p.m1)
  }
  ##########
  a.0=0.0001
  b.0=1
  a.1=0.5
  b.1=0.5
  a.2=rep(0.1,L)
  t_nub=15000
  burn=5000
  thin=1
  w.l=matrix(0,n,L)
  sigma.b1=rep(0,L)
  mu.b1=rep(0,L)
  w.lx=matrix(0,n,L)
  sigma.b1x=rep(0,L)
  mu.b1x=rep(0,L)
  fcn.1=rep(0,L)
  fco.1=rep(0,L)
  w.g=matrix(0,n,L)
  sigma.g1=rep(0,L)
  mu.g1=rep(0,L)
  p.a2=rep(0,p)
  w_e=matrix(0,n,L)
  b.e1=rep(0,L)
  w_c=matrix(0,n,L)
  p.l1=rep(0,L)
  p.l=rep(0,L)
  n.l=rep(0,L)
  #######
  beta.n=matrix(0,p,(t_nub-burn)/thin)
  b.n=matrix(0,L,(t_nub-burn)/thin)
  gama.n=matrix(0,p,(t_nub-burn)/thin)
  tb.g1s=matrix(0,p,(t_nub-burn)/thin)
  tb.g0s=matrix(0,p,(t_nub-burn)/thin)
  sigma.n=rep(0,(t_nub-burn)/thin)
  C.in=matrix(0,n,(t_nub-burn)*L/thin)
  w.n=matrix(0,L,(t_nub-burn)/thin)
  for(i in 1:t_nub){
    
    ###############beta.c
    for(b.j in 1:p){
      beta.sum=beta.c
      if(gama.c[b.j]==1){
        for(i.2 in 1:L){
          w.l[,i.2]=rep(1-tau[i.2],n)
          w.l[,i.2][y>b.c[i.2]+x%*%beta.c]=tau[i.2]
        }
        for(i.3 in 1:L){
          sigma.b1[i.3]=sum(2*C.im[,i.3]*w.l[,i.3]*(x[,b.j]^2)/sigma.c)
          mu.b1[i.3]=sum(2*C.im[,i.3]*w.l[,i.3]*x[,b.j]*(y-b.c[i.3]-x%*%beta.c+x[,b.j]*beta.c[b.j])/sigma.c)
        }
        sigma.b=(sum(sigma.b1)+1/delta.c[b.j])^(-1)
        mu.b=sigma.b*sum(mu.b1)
        beta.new=rnorm(1,mu.b,sqrt(sigma.b))
        beta.sum[b.j]=beta.new
        for(i.4 in 1:L){
          w.lx[,i.4]=rep(1-tau[i.4],n)
          w.lx[,i.4][y>b.c[i.4]+x%*%beta.sum]=tau[i.4]
        }
        for(i.5 in 1:L){
          sigma.b1x[i.5]=sum(2*C.im[,i.5]*w.lx[,i.5]*(x[,b.j]^2)/sigma.c)
          mu.b1x[i.5]=sum(2*C.im[,i.5]*w.lx[,i.5]*x[,b.j]*(y-b.c[i.5]-x%*%beta.sum+x[,b.j]*beta.new)/sigma.c)
        }
        sigma.bx=(sum(sigma.b1x)+1/delta.c[b.j])^(-1)
        mu.bx=sigma.bx*sum(mu.b1x)
        for(i.6 in 1:L){
          fcn.1[i.6]=sum(C.im[,i.6]*w.l[,i.6]*((y-b.c[i.6]-x%*%beta.sum)^2)/sigma.c)
          fco.1[i.6]=sum(C.im[,i.6]*w.lx[,i.6]*((y-b.c[i.6]-x%*%beta.c)^2)/sigma.c)
        }
        fcn=-sum(fcn.1)-beta.new^2/(2*delta.c[b.j])
        fco=-sum(fco.1)-beta.c[b.j]^2/(2*delta.c[b.j])
        fc=fcn-fco
        pro=-(beta.c[b.j]-mu.bx)^2/(2*sigma.bx)
        prn=-(beta.new-mu.b)^2/(2*sigma.b)
        prop=pro-prn
        d.a=(sigma.b/sigma.bx)^(0.5)*exp(fc+prop)
        if(runif(1)<d.a){
          beta.c[b.j]=beta.new
        }
      }
      if(gama.c[b.j]==0){
        beta.c[b.j]=0
      }
    }
    
    #######b.c
    for(i.b in 1:L){
      w_0=rep(1-tau[i.b],n)
      w_0[y>b.c[i.b]+x%*%beta.c]=tau[i.b]
      sigma.01=(sum(2*C.im[,i.b]*w_0/sigma.c))
      if(sigma.01>0){
        sigma.0=(sigma.01)^(-1)
        mu.0=sigma.0*sum(2*C.im[,i.b]*w_0*(y-x%*%beta.c)/sigma.c)
        b.new=rnorm(1,mu.0,sqrt(sigma.0))
        w_0x=rep(1-tau[i.b],n)
        w_0x[y>b.new+x%*%beta.c]=tau[i.b]
        sigma.0x=(sum(2*C.im[,i.b]*w_0x/sigma.c))^(-1)
        mu.0x=sigma.0x*sum(2*C.im[,i.b]*w_0x*(y-x%*%beta.c)/sigma.c)
        fcn.0=-sum(C.im[,i.b]*w_0*((y-b.new-x%*%beta.c)^2)/sigma.c)
        fco.0=-sum(C.im[,i.b]*w_0x*((y-b.c[i.b]-x%*%beta.c)^2)/sigma.c)
        fc.0=fcn.0-fco.0
        pro.0=-(b.c[i.b]-mu.0x)^2/(2*sigma.0x)
        prn.0=-(b.new-mu.0)^2/(2*sigma.0)
        prop.0=pro.0-prn.0
        d.0=(sigma.0/sigma.0x)^(0.5)*exp(fc.0+prop.0)
        if(runif(1)<d.0){
          b.c[i.b]=b.new
        }
      }
    }
    ###########gama.c
    for(a.j in 1:p){
      p.0=(p-sum(gama.c)+gama.c[a.j])
      p.1=(1+sum(gama.c)-gama.c[a.j])
      for(i.7 in 1:L){
        w.g[,i.7]=rep(1-tau[i.7],n)
        w.g[,i.7][y>b.c[i.7]+x%*%beta.c]=tau[i.7]
      }
      for(i.8 in 1:L){
        sigma.g1[i.8]=sum(2*C.im[,i.8]*w.g[,i.8]*(x[,a.j]^2)/sigma.c)
        mu.g1[i.8]=sum(2*C.im[,i.8]*w.g[,i.8]*x[,a.j]*(y-b.c[i.8]-x%*%beta.c+x[,a.j]*beta.c[a.j])/sigma.c)
      }
      h.a1=exp(-0.5*(sum(mu.g1))^2*((sum(sigma.g1)+1/delta.c[a.j])^(-1)))
      h.a2=(1+delta.c[a.j]*sum(sigma.g1))^(0.5)
      h.a=h.a1*h.a2*p.0/p.1
      p.a2[a.j]=1/(1+h.a)
      gama.c[a.j]=rbinom(1,1,p.a2[a.j])
    }
    
    ##########sigma.c
    for(i.e in 1:L){
      w_e[,i.e]=rep(1-tau[i.e],n)
      w_e[,i.e][y>b.c[i.e]+x%*%beta.c]=tau[i.e]
    }
    a.e1=rowSums(C.im)
    a.e=a.0+sum(a.e1)/2
    for(i.e1 in 1:L){
      b.e1[i.e1]=sum(C.im[,i.e1]*w_e[,i.e1]*(y-b.c[i.e1]-x%*%beta.c)^2)
    }
    b.e=b.0+sum(b.e1)
    sigma.c=rinvgamma(1,a.e,b.e)
    
    ########delta.c
    for(i.d in 1:p){
      if(gama.c[i.d]==1){
        delta.c[i.d]=rinvgamma(1,a.1+0.5,b.1+(beta.c[i.d]^2)/2)
      }
      if(gama.c[i.d]==0){
        delta.c[i.d]=rinvgamma(1,a.1,b.1)
      }
    }
    
    #######w.k
    for(w.1 in 1:L){
      n.l[w.1]=sum(C.im[,w.1])
    }
    n.lx=n.l+a.2
    w.k=rdirichlet(1,n.lx)
    
    ###########C.im
    for(i.c2 in 1:L){
      w_c[,i.c2]=rep(1-tau[i.c2],n)
      w_c[,i.c2][y>b.c[i.c2]+x%*%beta.c]=tau[i.c2]
    }
    for(i.c in 1:n){
      for(i.c1 in 1:L){
        p.l1[i.c1]=w.k[i.c1]*exp(-w_c[i.c,i.c1]*(y[i.c]-b.c[i.c1]-x[i.c,]%*%beta.c)^2/sigma.c)
      }
      for(i.9 in 1:L){
        p.l[i.9]=p.l1[i.9]/sum(p.l1)
      }
      C.im[i.c,]=rmultinom(1,1,p.l)
    }
    
    ################
    tb.g1=rep(0,p)
    for(i.pb in 1:p){
      if(gama.c[i.pb]==1){
        tb.g1[i.pb]=tb.g1[i.pb]+1
      }
      
    }
    ##################
    tb.g0=rep(0,p)
    for(i.pb1 in 1:p){
      if(gama.c[i.pb1]==0){
        tb.g0[i.pb1]=tb.g0[i.pb1]+1
      }
      
    }
    
    if((i-burn)/thin%in%1:(t_nub-burn)/thin){
      ib=(1+((i-burn)/thin-1)*L):((i-burn)*L/thin)
      tb.g1s[,(i-burn)/thin]=tb.g1
      tb.g0s[,(i-burn)/thin]=tb.g0
      gama.n[,(i-burn)/thin]=gama.c
      b.n[,(i-burn)/thin]=b.c
      beta.n[,(i-burn)/thin]=beta.c
      sigma.n[(i-burn)/thin]=sigma.c
      C.in[,ib]=C.im
      w.n[,(i-burn)/thin]=w.k
    }
    
  }

  beta.m=apply(beta.n,1,mean)
  y.gu=x%*%beta.m
  sigma.m=mean(sigma.n)
  burn=5000
  mcmc=15000
  DIC1=dicli(y, x, beta.n, beta.m, b.n,
             burn, mcmc,sigma.n,C.in,w.n) 
  MSE=mean((y.gu-y)^2)
  RMSE=sqrt(MSE)
  AD=mean(abs(y-y.gu))
  MAD.ii=rep(0,10000)
  for(i in 1:10000){
    y.gui=x%*%beta.n[,i]
    MAD.ii[i]=mean(abs(y-y.gui))
  }
  mean(MAD.ii)
  
 sd((y.gu-y)^2)
  tb.g0m=apply(tb.g0s,1,mean)
  tb.g1m=apply(tb.g1s,1,mean)
  for(i.tb in 1:p){
    tt.g[i.tb]=max(tb.g0m[i.tb],tb.g1m[i.tb])
    if(tt.g[i.tb]==tb.g0m[i.tb]){
      z.v=z.v+1
    }
    if(tt.g[i.tb]==tb.g1m[i.tb]){
      l.v=l.v+1
    }
  }
  
  RMSE.c[i.x]=RMSE
  beta.mm[,i.x]=beta.m
}
})


