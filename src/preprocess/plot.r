a <- read.table("groundTruth_G1")
a <- a[(a$V3+a$V4)>=40,]
n=50;
ta <- a[(a$V3+a$V4)==n & a$V8!=2 & a$V10!=a$V11 & a$V8!=a$V9,]; histc <- hist(ta$V3, col="blue", breaks=seq(-0.5,n+0.5,by=1))
print(mean(ta$V3))
print(sd(ta$V3))
tb <- a[(a$V3+a$V4)==n & a$V8!=2 & a$V10==a$V11 & a$V8!=a$V9 & a$V10==1,]; histb <- hist(tb$V3, col="blue", breaks=seq(-0.5,n+0.5,by=1))
print(mean(tb$V3))
print(sd(tb$V3))
tc <- a[(a$V3+a$V4)==n & a$V8!=2 & a$V10==a$V11 & a$V8!=a$V9 & a$V10==0,]; hista <- hist(tc$V3, col="blue", breaks=seq(-0.5,n+0.5,by=1))
print(mean(tc$V3))
print(sd(tc$V3))
td <- a[(a$V3+a$V4)==n & a$V8!=2 & a$V10!=a$V11 & a$V8==a$V9,]; histd <- hist(td$V3, col="blue", breaks=seq(-0.5,n+0.5,by=1))
te <- a[(a$V3+a$V4)==n & a$V8!=2 & a$V10==a$V11 & a$V8==a$V9,]; histe <- hist(te$V3, col="blue", breaks=seq(-0.5,n+0.5,by=1))
tf <- a[(a$V3+a$V4)==n & a$V8==2,]; histf <- hist(tf$V3, col="blue", breaks=seq(-0.5,n+0.5,by=1))

#x <- a[a$V3/(a$V3+a$V4)>0.2 & a$V3/(a$V3+a$V4)<0.8,]
#hist(x$V3/(x$V3+x$V4), col="blue", breaks=seq(-0.01,1.01,by=0.02))


#ta <- a[a$V8!=2 & a$V10!=a$V11 & a$V8!=a$V9,]; histc <- hist(ta$V3/(ta$V3+ta$V4), col="blue", breaks=seq(-0.01,1.01,by=0.02))
#tb <- a[a$V8!=2 & a$V10==a$V11 & a$V8!=a$V9 & a$V10==1,]; histb <- hist(tb$V3/(tb$V3+tb$V4), col="blue", breaks=seq(-0.01,1.01,by=0.02))
#tc <- a[a$V8!=2 & a$V10==a$V11 & a$V8!=a$V9 & a$V10==0,]; hista <- hist(tc$V3/(tc$V3+tc$V4), col="blue", breaks=seq(-0.01,1.01,by=0.02))
#td <- a[a$V8!=2 & a$V10!=a$V11 & a$V8==a$V9,]; histd <- hist(td$V3/(td$V3+td$V4), col="blue", breaks=seq(-0.01,1.01,by=0.02))
#te <- a[a$V8!=2 & a$V10==a$V11 & a$V8==a$V9,]; histe <- hist(te$V3/(te$V3+te$V4), col="blue", breaks=seq(-0.01,1.01,by=0.02))
#tf <- a[a$V8==2,]; histf <- hist(tf$V3/(tf$V3+tf$V4), col="blue", breaks=seq(-0.01,1.01,by=0.02))

print(dim(ta)[1])
ca=dim(ta)[1]
cb=dim(tb)[1]
cc=dim(tc)[1]
cd=dim(td)[1]
total=ca+cb+cc+cd
print(ca/total)
print(cb/total)
print(cc/total)
print(cd/total)
print(dim(td[td$V8==0,])[1]/total)
print(dim(td[td$V8==1,])[1]/total)
plot(histd,col=rgb(1,1,0,0.5))
lines(histb,col=rgb(0,1,0,0.5))
lines(histc,col=rgb(0,0,1,0.5))
lines(hista,col=rgb(1,0,0,0.5))
lines(histe,col=rgb(0,0,0,0.5))
