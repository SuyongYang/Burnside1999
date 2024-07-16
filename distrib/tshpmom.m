function [facv,facr]=tshpmom(sigma,drules,ncorr) ;

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
hbig=[ eye(ns+nex) ; h ] ;

base=181;
maxf=101;

jp=(1:1:maxf)' ;
jm=(maxf:-1:1)' ;
hpap=-(0.894.^jp).*(0.0561*cos(jp*0.112)+0.0558*sin(jp*0.112)) ;
hpam=-(0.894.^jm).*(0.0561*cos(jm*0.112)+0.0558*sin(jm*0.112)) ;
hpa=[ hpam ; 1-(0.0561*cos(0)+0.0558*sin(0)) ; hpap ] ;

gammf=kron(zeros(1,ncorr+1),zeros(ns+nex,ns+nex)) ;
facv=zeros(ny,(ncorr+1)*ny) ;

for k=0:ncorr
  for j=0:base
    if j==0
      gammj=gamm0 ;
      gammf(:,k*(ns+nex)+1:(k+1)*(ns+nex))=gammf(:,k*(ns+nex)+1:(k+1)*(ns+nex))+gammj*(hpa(k+1:2*maxf+1,1)'*hpa(1:2*maxf+1-k,1)) ;
    else
      gammj=m*gammj ;
      if j<=k
        gammf(:,k*(ns+nex)+1:(k+1)*(ns+nex))=gammf(:,k*(ns+nex)+1:(k+1)*(ns+nex))+gammj*(hpa(k+1-j:2*maxf+1,1)'*hpa(1:2*maxf+1-k+j,1))+gammj'*(hpa(k+j+1:2*maxf+1,1)'*hpa(1:2*maxf+1-k-j,1)) ;
      else
        gammf(:,k*(ns+nex)+1:(k+1)*(ns+nex))=gammf(:,k*(ns+nex)+1:(k+1)*(ns+nex))+gammj*(hpa(1:2*maxf+1+k-j,1)'*hpa(1-k+j:2*maxf+1,1))+gammj'*(hpa(k+j+1:2*maxf+1,1)'*hpa(1:2*maxf+1-k-j,1)) ;
      end
    end
  end
  facv(:,k*ny+1:(k+1)*ny)=hbig*gammf(:,k*(ns+nex)+1:(k+1)*(ns+nex))*(hbig') ;
end

sd=sqrt(diag(facv(1:ny,1:ny))) ;
tcorr=kron(sd,sd') ;
facr=facv./kron(ones(1,ncorr+1),tcorr) ;
