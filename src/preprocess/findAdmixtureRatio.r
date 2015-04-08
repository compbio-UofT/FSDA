args<-commandArgs(TRUE)
s <- read.table(args[1])

s <- s[ s$V3/(s$V3+s$V4)<0.2,]
p<- data.frame(r=numeric(0),l=numeric(0))
#hist(s$V2/(s$V2+s$V3),col="blue")
#q()
step=0.0005

for (r in seq(0.01, 0.1, by=step) ) {
	#print(r)
	ll=0;
	for (i in seq(10,61))
	{
		k=100
		x=s[s$V3+s$V4==i,] 
		h <- hist(x$V3, breaks=seq(0,k,by=1), plot=FALSE);
		k=1
		while (h$counts[k]>0)
		{
			k=k+1;
		}
		y=dbinom(seq(1,k),i,r)
		y=y/(sum(y))
		y=log(y)
		ll=ll+sum(y*h$counts)
	}
	p=rbind(p,data.frame(r=2*r,l=ll))
}
p$l=p$l-max(p$l)
p$l=exp(p$l)
p$l=p$l/sum(p$l)
p$l=p$l/(2*step)
write(p$r,file=args[2], ncolumns=length(p$l))
write(p$l,file=args[2], ncolumns=length(p$l), append=TRUE)
#plot(p$r,p$l,type="l",xlim=c(0.05,0.16), col="blue")

