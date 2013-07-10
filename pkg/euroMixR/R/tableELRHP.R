tableELRHP=function(L=4){
res=NULL
p=rep(1,L)
p=p/sum(p)
q=q012(p=p)
probHD=qkappa(kappa=c(1,0.0,0.0),q=q)
k0=c(1,0,0,0.25,0.5,0.75,9/16,17/32)
k1=c(0,1,0,0.5,0.5,0.25,6/16,14/32)
k2=c(0,0,1,0.25,0,0,1/16,1/32)
phi=(2*k2+k1)/4
kappaMat=cbind(k0,k1,k2)
for (i in 1:dim(kappaMat)[1]){
k0=as.double(kappaMat[i,1]);k1=as.double(kappaMat[i,2]);k2=as.double(kappaMat[i,3])
exact=k0^2+2*k1*k0+k1^2*(L+3)/4+2*k2*k0+k2*k1*(L+1)+k2^2*L*(L+1)/2
b=(k1^2+4*k1*k2+2*k2^2)/4
c=k2^2/2
a=1-b-c
exact2=a+b*L+c*L^2
probHP=qkappa(kappa=c(k0,k1,k2),q=q)
LR.num=sum((probHP)^2/probHD)
ledd=c(k0=k0,k1=k1,k2=k2,formel2=exact2,num=LR.num)
res=rbind(res,ledd)
}
nn=c("Unrelated","Parent-offspring","Monozygous-twin","Full-Sib","Half-sib etc","First cousin","Double first-cousin","Quadruple...")
rownames(res)=nn
res
}
