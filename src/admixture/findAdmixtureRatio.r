library("VGAM")
args<-commandArgs(TRUE)

s <- read.table(args[1])
n=40
s <- s[ (s$V3+s$V4)>=n & s$V3/(s$V3+s$V4)<0.2,]
print(mean(s$V3/(s$V3+s$V4))*2)
print(dim(s))


p<- data.frame(r=numeric(0),l=numeric(0))
step=0.0005
covs=seq(40,100)
hists=c()
dropout=TRUE
if (dropout)
	s <- s[s$V3>0,]
print(mean(s$V3/(s$V3+s$V4))*2)
for (i in covs)
{
	x=s[s$V3+s$V4==i,] 
	h <- hist(x$V3, breaks=seq(-1,i,by=1), plot=FALSE);
	hists=c(hists,list(h))
}



for (r in seq(0.01, 0.1, by=step) ) {
	#print(r)
	ll=0;
	ll_list=c()
	for (l in seq(0,0.02,by=0.002))
	{
		tmpll=0
		for (i in covs)
		{
			y=dbetabinom(seq(0,i), i, r, rho = l, log = FALSE)
			if (dropout)
				y[1]=0 ##
			y=y/(sum(y))
			y=log(y)
			if (dropout)
				y[1]=0 ##
			tmpll=tmpll+sum(y*hists[[i-covs[1]+1]]$counts)
		}
		ll_list=c(ll_list,tmpll)
	}
	mx=max(ll_list)
	ll=mx+log(sum(exp(ll_list-mx)))
	p=rbind(p,data.frame(r=2*r,l=ll))
}

p$l=p$l-max(p$l)
p$l=exp(p$l)
p$l=p$l/sum(p$l)
p$l=p$l/(2*step)
bestR=r[1]
bestL=l[1]
for (i in seq(1,length(p$l)) ) {
	if (p$l[i]>0.1)
	{
#		print (p$r[i])
#		print (p$l[i])
		if (bestL<p$l[i])
		{
			bestL=p$l[i]
			bestR=p$r[i]
		}
	}
}
print("_____________")
print (bestR)
bestL=0
bestRHO=0

for (l in seq(0,0.02,by=0.002))
	{
		tmpll=0
		for (i in covs)
		{
			y=dbetabinom(seq(0,i), i, bestR/2, rho = l, log = FALSE)
			if (dropout)
				y[1]=0 ##
			y=y/(sum(y))
			y=log(y)
			if (dropout)
				y[1]=0 ##
			tmpll=tmpll+sum(y*hists[[i-covs[1]+1]]$counts)
			#totalsum+=tmpll
		}
		if (l==0 | tmpll>bestL)
		{
			bestL=tmpll
			bestRHO=l
		}
	}
print (bestRHO)

plot(p$r,p$l,type="l",xlim=c(0.05,0.16), col="blue")


#s <- s[ (s$V3+s$V4)==50,]
h <- hist(s$V3,col="blue", breaks=seq(-0.5,50.5,by=1), freq=FALSE,plot=FALSE)
x=seq(0,n)

y=dbetabinom(x, n, bestR/2, rho = bestRHO, log = FALSE)
plot(h,freq=FALSE, col="blue",ylim=c(0,0.35))
points(x,y,  col="red")
q()

