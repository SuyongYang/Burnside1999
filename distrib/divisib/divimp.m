function irf=divimp(b,nimp)

beta_=1.03^(-.25) ;
sigma=b(4,1) ; sigma=sigma*sigma ;

drules=divsolve(beta_,b(1,1),b(8,1),b(6,1),b(7,1),b(3,1)) ;
[nex,nex2]=size(sigma) ;
[ny,ns]=size(drules) ;
ns=ns-nex ;
m=drules(1:nex+ns,:) ;
h=drules(nex+ns+1:ny,:) ;
irf=zeros(nimp,9) ;

for k=1:nimp
  if k==1
    irfs=eye(ns+nex) ;
  else
    irfs=m*irfs ;
  end
  irff=h*irfs ;
  irf(k,1)=irfs(ns+1,ns+1) ;
  irf(k,2)=irfs(1,ns+1) ;
  irf(k,3)=irff(2,ns+1) ;
  irf(k,4)=irff(3,ns+1) ;
  irf(k,5)=irff(4,ns+1) ;
  irf(k,6)=irff(5,ns+1) ;
  irf(k,7)=irff(6,ns+1) ;
  irf(k,8)=irff(7,ns+1) ;
  irf(k,9)=irff(8,ns+1) ;
end

t=(1:1:nimp)' ;
zers=zeros(nimp,1) ;

figure(1)

subplot(221)
plot(t,[ irf(:,1) zers],'k-' )
title('A')
set(gca,'XLim',[0 nimp+1]) ;

subplot(222)
plot(t,[ irf(:,4) zers ],'k-',t,irf(:,6),'k:')
title('H and Y/N')
set(gca,'XLim',[0 nimp+1]) ;

subplot(223)
plot(t,[ irf(:,3) zers ],'k-',t,irf(:,7),'k:')
title('C and I')
set(gca,'XLim',[0 nimp+1]) ;

subplot(224)
plot(t,[ irf(:,5) zers ],'k-',t,irf(:,2),'k:')
title('Y and K')
set(gca,'XLim',[0 nimp+1]) ;
