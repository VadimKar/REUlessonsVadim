

library(MASS)

mf=function(xcor=c(L=0.5,H=0.5),k=12,xms=c(L=0.25,H=0.75),sgn=1,sdx=0.25){
	# https://math.stackexchange.com/questions/2496649/derivative-and-second-derivative-of-scaled-logistic
	#Equations of the sigmoid and its first 2 derivatives
	F=function(x) 1/(1+exp(-k*(x-0.5)))
	dF=function(x) k*exp(-k*(x-0.5))/((1+exp(-k*(x-0.5)))^2)
	d2F=function(x) k*k*(exp(k*(x-0.5))-1)*exp(k*(x-0.5))/((exp(k*(x-0.5))+1)^3)
	
	SF=function(rho) matrix(c(1,rho*sqrt(pi/3),rho*sqrt(pi/3),1), 2,2, byrow=TRUE)
	set.seed(1); XL=mvrnorm(1e4,c(0,0),Sigma=SF(xcor[1]))
	set.seed(1); XH=mvrnorm(1e4,c(0,0),Sigma=SF(xcor[2]))
	YL=cor(F(xms[1]+sdx*XL*sgn))[1,2]
	YH=cor(F(xms[2]+sdx*XH*sgn))[1,2]
	sim=c(YL-YH,YL,YH)
	
	xcov=xcor*sdx^2
	ycori=function(ind) (xcor[ind] + 0.5*(d2F(xms[ind])^2)*xcor[ind]*xcov[ind]*(dF(xms[ind])^-2)) / (1+0.5*(d2F(xms[ind])^2)*sdx*sdx*(dF(xms[ind])^-2))
	an=c(ycori(1)-ycori(2),ycori(1),ycori(2))
	xs=seq(0,1,len=1e2); plot(xs,F(xs)); abline(v=xms,lwd=2,col=c(2,4)); abline(v=xms-sdx,col=c(2,4),lty=2); abline(v=xms+sdx,col=c(2,4),lty=2);
	round(rbind(sim,an),2)
}
#Some examples
mf(sdx=0.15,k=18,xms=c(.3,.6),xcor=c(.3,.8))
mf(sdx=0.15,k=25,xms=c(.2,.5),xcor=rep(.9,2),sgn=-1)
mf(sdx=0.1,k=25,xms=c(.4,.7),xcor=rep(.9,2),sgn=1)


