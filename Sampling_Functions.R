library(truncnorm)
library(invgamma)

beta.samp=function(x,sig,R.inv,b){
	v=solve(t(x)%*%x+(1/sig)*R.inv)
	mu=v%*%t(x)%*%(Y-b)
	res=as.vector(1,mu,v)
	return(res)
}

sig.samp=function(y,beta,x,b,as,bs){
	n=length(y)
	a=(n/2)+as
	b=(1/2)*sum((y-x%*%beta-b)^2)+bs
	res=1/rgamma(1,shape=a,rate=b)
	return(res)
}

tau.samp=function(b,rho,D,W,at,bt){
	n=dim(D)[1]
	a=(n/2)+at
	b=(1/2)*t(b)%*%(D-rho*W)%*%b+bt
	res=rinvgamma(1,shape=a,rate=b)
	return(res)
}

b.samp=function(y,x,sig,tau,rho,D,W){
	n=dim(D)[1]
	Q = (1/sig)*diag(n) + (1/tau)*(D-rho*W)
	b = (y-x%*%beta)*(1/sig)
	L = t(chol(Q))
	v = forwardsolve(L,b)
	mu = backsolve(t(L),v)
	z = rnorm(n,0,1)
	y = backsolve(t(L),z)
	res = y + mu
	return(res)
}

rho.samp=function(rho.i,b,tau,D,W,lb,ub){
  r.n=rtruncnorm(1,a=lb,b=ub,mean=rho.i,sd=0.05)
  r.n.p=(-1/2)*log(det(D-r.n*W))+(1/2)*r.n*(1/tau)*t(b)%*%W%*%b
  r.n.ji=log(dtruncnorm(rho.i,a=lb,b=ub,mean=r.new,sd=0.05))
  r.i.p=(-1/2)*log(det(D-rho.i*W))+(1/2)*rho.i*(1/tau)*t(b)%*%W%*%b
  r.i.jn=log(dtruncnorm(r.n,a=lb,b=ub,mean=rho.i,sd=0.05))
  ratio=exp(r.n.p+r.n.ji-r.i.p-r.i.jn)
  ratio=min(1,ratio)
  q=rbinom(1,1,ratio)
  res=q*r.n+(1-q)*rho.i
  return(list(res,q))
}

  
  
	
	