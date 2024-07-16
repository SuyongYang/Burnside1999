function mh=bersolve(theta_,alpha_,deltak,gammax,sg,rho)

global beta_ capt xi ssf 

nc=2 ;
ns=2 ;
ncs=1 ;
nex=2 ;
nf=4 ;

gammax=exp(gammax) ;
kyratio=(1-alpha_)/(gammax/beta_-(1-deltak)) ;
si=(gammax-(1-deltak))*kyratio ;
sc=1-si-sg ;

etaa=(gammax-beta_*(1-deltak))/gammax;
etak=-alpha_*etaa ;
etan=alpha_*etaa ;

cons0=theta_*log(capt/(capt-xi-ssf)) ;
cons1=ssf/(capt-xi-ssf) ;
cons2=theta_*cons1/cons0 ;

mcc = [ -1    0
         0 1-alpha_+cons1 ] ;
mcs = [ 0       0        1
       1-alpha_ alpha_-1 1 ] ;
mce = [ 0 0
        1 0 ] ;
mss0 = [ -etak           -etan    -1
         alpha_-1       1-alpha_  -1
        -gammax*kyratio    0       0 ]  ;
mss1 = [ 0                             0     1
         0                             0     0
        1-alpha_+(1-deltak)*kyratio  alpha_  0 ] ;
msc0 = [ 0     etan
         0 alpha_-cons2
         0       0      ] ;
msc1 = [  0     0
          0     0
         sc  -alpha_ ] ;
mse0 = [ etaa 0
          1   0
          0   0  ] ;
mse1 = [ 0   0
         0   0
        -1  sg ] ;

fc = [ 0     alpha_
       0     alpha_
     -sc/si alpha_/si
       0     alpha_    ] ;
fx = [ 1-alpha_         alpha_
       1-alpha_        alpha_-1
       (1-alpha_)/si  alpha_/si
          0                0    ] ;
fe = [ 1      0
       1      0
      1/si  -sg/si
       1      0  ] ;

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
solc=[ cxl(:,1:ns) ce ] + cxl(:,ns+1:ns+ncs)*soll ;

solf=[ fx fe ]+fc*solc ;

m = [ solx ; [ zeros(nex,ns) rho ] ] ;
h = [ soll ; solc ; solf ] ;
mh=[ m ; h ] ;
