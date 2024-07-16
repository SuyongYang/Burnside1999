function mh=cesolve(beta_,theta_,alpha_,gamm,deltak,g,rho) ;

nc=2 ;
ns=1 ;
ncs=1 ;
nex=2 ;
nf=5 ;

kyratio=(1-alpha_)/(1/beta_-(1-deltak)*exp(-gamm)) ;
knratio=((1/beta_-(1-deltak)*exp(-gamm))/((1-alpha_)*exp((alpha_-1)*gamm)))^(-1/alpha_) ;
ssc=alpha_*(knratio^(1-alpha_))*exp((alpha_-1)*gamm)/theta_ ;
ssy=(ssc+g)/(1-(1-(1-deltak)*exp(-gamm))*kyratio) ;

iyratio=(1-(1-deltak)*exp(-gamm))*kyratio ;
cyratio=ssc/ssy ;

mu=beta_*(1-alpha_)/kyratio ;
muk=-alpha_*mu ;
mun=alpha_*mu ;

mcc = [ -1 0
         0 1-alpha_ ] ;
mcs = [ 0        1
        1-alpha_ 1 ] ;
mce = [ 0       0
       alpha_-1 0 ] ;
mss0 = [ muk    1
       -kyratio 0 ] ;
mss1 = [ 0                                    -1
       1-alpha_+(1-deltak)*exp(-gamm)*kyratio  0 ] ;
msc0 = [ 0  -mun
         0  0 ] ;
msc1 = [ 0         0
       cyratio  -alpha_ ] ;
mse0 = [ 1-alpha_*mu 0
         0           0 ] ;
mse1 = [ 0                                         0
       (1-alpha_)+(1-deltak)*exp(-gamm)*kyratio  g/ssy ] ;

fc = [ 0 alpha_
       0 alpha_-1
     -cyratio/iyratio alpha_/iyratio
       0 alpha_*mu
      -1  0 ] ;
fx = [ 1-alpha_
       1-alpha_
     (1-alpha_)/iyratio
      -alpha_*mu
         0 ] ;
fe = [ alpha_-1  0
       alpha_-1  0
     (alpha_-1)/iyratio  -g/(ssy*iyratio)
      alpha_*mu  0
         0       0 ] ;

w = -inv(mss0 - msc0*inv(mcc)*mcs)*(mss1 - msc1*inv(mcc)*mcs);
r =  inv(mss0 - msc0*inv(mcc)*mcs)*(mse0 + msc0*inv(mcc)*mce);
q =  inv(mss0 - msc0*inv(mcc)*mcs)*(mse1 + msc1*inv(mcc)*mce);

[pr,lambr]=eig(w);
alamb=abs(diag(lambr)) ;
[lambs,lambz]=sort(alamb) ;

lambda=lambr(lambz,lambz) ;
p=pr(:,lambz) ;

lamb1=lambda(1:ns,1:ns) ;
lamb2=lambda(ns+1:ns+ncs,ns+1:ns+ncs) ;

p11=p(1:ns,1:ns) ;
p12=p(1:ns,ns+1:ns+ncs) ;
p21=p(ns+1:ns+ncs,1:ns) ;
p22=p(ns+1:ns+ncs,ns+1:ns+ncs) ;

ps=inv(p) ;
ps11=ps(1:ns,1:ns) ;
ps12=ps(1:ns,ns+1:ns+ncs) ;
ps21=ps(ns+1:ns+ncs,1:ns) ;
ps22=ps(ns+1:ns+ncs,ns+1:ns+ncs) ;

rxe=r(1:ns,1:nex) ;
rle=r(ns+1:ns+ncs,1:nex) ;
qxe=q(1:ns,1:nex) ;
qle=q(ns+1:ns+ncs,1:nex) ;

phi0=ps21*rxe+ps22*rle ;
phi1=ps21*qxe+ps22*qle ;

psi=zeros(ncs,nex) ;

for i=1:ncs
 psi(i,:)=-(phi0(i,:)*rho+phi1(i,:))*inv(eye(nex)-rho/lamb2(i,i))/lamb2(i,i) ;
end

xx=p11*lamb1*inv(p11) ;
xe=(p11*lamb1*ps12+p12*lamb2*ps22)*inv(ps22)*psi+rxe*rho+qxe ;
solx=[ xx xe ] ;

lx=-inv(ps22)*ps21 ;
lex=inv(ps22)*psi ;
soll=[lx lex] ;

cxl=inv(mcc)*mcs ;
ce=inv(mcc)*mce ;
solc=[ cxl*[ eye(ns) ; lx ] cxl*[ zeros(ns,nex) ; lex ]+ce ] ;

% update the interest rate equation

fx(5,:)=solc(1,1:ns)*xx ;
fe(5,:)=solc(1,1:ns)*xe+solc(1,ns+1:ns+nex)*rho ;

solf=[ fx fe ]+fc*solc ;

m = [ solx ; [ zeros(nex,ns) rho ] ] ;
h = [ soll ; solc ; solf ] ;
mh=[ m ; h ] ;

