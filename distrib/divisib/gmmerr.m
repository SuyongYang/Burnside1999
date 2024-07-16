function u=gmmerr(b)

global xdata ;

dep=1-(xdata(2:115,5)-xdata(1:114,2))./xdata(1:114,5) ;

t=(1:1:113)' ;
t0=(0:1:114)' ;

beta_=1.03^(-.25) ;
loga=log(xdata(1:115,3))-log(xdata(1:115,5))*(1-b(8,1))-log(xdata(1:115,4))*b(8,1)-t0*b(6,1)*b(8,1) ;

u1=(ones(113,1)*b(1,1))./(ones(113,1)-xdata(2:114,4))-b(8,1)*xdata(2:114,3)./(xdata(2:114,1).*xdata(2:114,4)) ;
u2=loga(2:114,1)-b(2,1)-loga(1:113,1)*b(3,1) ;
u3=u2.*loga(1:113,1) ;
u4=u2.^2-b(4,1)^2 ;
u5=log(xdata(2:114,3))-b(5,1)-t*b(6,1) ;
u6=u5.*t/113 ;
u7=dep(2:114,1)-b(7,1) ;
u8=1-beta_*(xdata(2:114,1)./xdata(3:115,1)).*((1-b(8,1))*xdata(3:115,3)./xdata(3:115,5)+1-b(7,1)) ;
u9=xdata(2:114,8).^2-b(9,1)^2 ;
u10=xdata(2:114,6).^2-b(10,1)^2*xdata(2:114,8).^2 ;
u11=xdata(2:114,7).^2-b(11,1)^2*xdata(2:114,8).^2 ;
u12=xdata(2:114,9).^2-b(12,1)^2*xdata(2:114,8).^2 ;
u13=xdata(2:114,9).^2-b(13,1)^2*xdata(2:114,10).^2 ;

u=[ u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 ] ;
