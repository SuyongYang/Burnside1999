function [facv,facr]=dshpmom(sigma,drules,indgro,ncorr) ;

[nex,nex2]=size(sigma) ;
[ny,ns]=size(drules) ;

ns=ns-nex ;

m=drules(1:nex+ns,:) ;
h=drules(nex+ns+1:ny,:) ;
sigh=zeros(ns+nex,ns+nex) ;
sigh(ns+1:ns+nex,ns+1:ns+nex)=sigma ;

[vr,dr]=eig(m) ;
dr=diag(dr) ;

sigt=inv(vr)*sigh*inv(vr') ;
gammt0=(ones(ns+nex,ns+nex)./(ones(ns+nex,ns+nex)-kron(dr,dr'))).*sigt ;
gamm0=vr*gammt0*(vr') ;

a0=m(1:ns,1:ns) ; a1=m(1:ns,ns+1:ns+nex) ;
pi0=m(ns+1:ns+nex,ns+1:ns+nex) ;
pibar=[ [pi0 ; eye(nex) ] zeros(2*nex,nex) ] ;
thx=[ ones(ns,1) zeros(ns,nex-1) ] ;
mbar=[ [ a0  a1+thx  (-(a1+a0*thx)) ] ; [ zeros(2*nex,ns) pibar ] ] ;
b0=[ [ eye(ns) zeros(ns,nex) ] ; [ zeros(nex,ns) eye(nex) ] ; zeros(nex,ns+nex) ] ;
b1=[ [ -eye(ns) thx ] ; zeros(nex,ns+nex) ; [ zeros(nex,ns) eye(nex) ] ] ;

gamm0=b0*gamm0*(b0')+b1*gamm0*(b1')+b0*m*gamm0*(b1')+(b0*m*gamm0*(b1'))' ;

nf=ny-nex-ns ;
h1=drules(nex+ns+1:nex+ns+nf,1:ns) ;
h2=drules(nex+ns+1:nex+ns+nf,ns+1:ns+nex) ;
thf=[ indgro zeros(nf,nex-1) ] ;
hbar=[ h1 h2+thf (-(h2+h1*thx)) ] ;
hbig=[ eye(ns) zeros(ns,2*nex) ;
 zeros(1,ns) 1 zeros(1,2*nex-1) ;
 zeros(nex-1,ns) ones(nex-1,1) eye(nex-1) zeros(nex-1,1) (-eye(nex-1)) ;
 hbar ] ;

base=401;
maxf=251;

jp=(1:1:maxf)' ;
jm=(maxf:-1:1)' ;
hpap=-(0.894.^jp).*(0.0561*cos(jp*0.112)+0.0558*sin(jp*0.112)) ;
hpam=-(0.894.^jm).*(0.0561*cos(jm*0.112)+0.0558*sin(jm*0.112)) ;
hpa=[ hpam ; 1-(0.0561*cos(0)+0.0558*sin(0)) ; hpap ] ;
hpa=cumsum(hpa) ;

facv=zeros(ny,(ncorr+1)*ny) ;

for k=0:ncorr
  gammf=zeros(ns+2*nex,ns+2*nex) ;
  for j=0:base
    if j==0
      gammj=gamm0 ;
      gammf=gammf+gammj*(hpa(k+1:2*maxf+1,1)'*hpa(1:2*maxf+1-k,1)) ;
    else
      gammj=mbar*gammj ;
      if j<=k
        gammf=gammf+gammj*(hpa(k+1-j:2*maxf+1,1)'*hpa(1:2*maxf+1-k+j,1))+gammj'*(hpa(k+j+1:2*maxf+1,1)'*hpa(1:2*maxf+1-k-j,1)) ;
      else
        gammf=gammf+gammj*(hpa(1:2*maxf+1+k-j,1)'*hpa(1-k+j:2*maxf+1,1))+gammj'*(hpa(k+j+1:2*maxf+1,1)'*hpa(1:2*maxf+1-k-j,1)) ;
      end
    end
  end
  facv(:,k*ny+1:(k+1)*ny)=hbig*gammf*(hbig') ;
end

sd=sqrt(diag(facv(1:ny,1:ny))) ;
tcorr=kron(sd,sd') ;
facr=facv./kron(ones(1,ncorr+1),tcorr) ;
