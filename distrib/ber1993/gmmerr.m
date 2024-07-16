function u=gmmerr(b)

global xdata beta_ capt xi ssf 

sg=exp(b(6,1)-b(4,1)) ;
rho=[ b(11,1) 0 ; 0 b(8,1) ] ;
ssn=getssn(b(1,1),b(2,1),b(3,1),b(5,1),sg) ;

dep=1-(xdata(2:115,6)-xdata(1:114,2))./xdata(1:114,6) ;
weff=(xdata(1:115,4)./(xdata(1:115,5).*xdata(1:115,1)*b(1,1)+b(2,1)*ssf*xdata(1:115,4)))*b(2,1)*(capt-xi) ;

t=(1:1:113)' ;
t0=(0:1:114)' ;

loga=log(xdata(1:115,4))-log(xdata(1:115,6))*(1-b(2,1))-log(xdata(1:115,5))*b(2,1)-log(weff)*b(2,1)-t0*b(5,1)*b(2,1) ;

u1=log(xdata(2:114,5))-log(ssf)-log(ssn) ;
u2=1-beta_*(xdata(2:114,1)./xdata(3:115,1)).*((1-b(2,1))*xdata(3:115,4)./xdata(3:115,6)+1-b(3,1)) ;
u3=dep(2:114,1)-b(3,1) ;
u4=log(xdata(2:114,4))-b(4,1)-t*b(5,1) ;
u5=u4.*t/113 ;
u6=log(xdata(2:114,3))-b(6,1)-t*b(7,1) ;
u7=u6.*t/113 ;
u6l=log(xdata(1:113,3))-b(6,1)-(t-1)*b(7,1) ;
u8=(u6-u6l*b(8,1)).*u6l ;
u9=(u6-u6l*b(8,1)).^2-b(9,1)^2 ;
u10=loga(2:114,1)-b(10,1)-loga(1:113,1)*b(11,1) ;
u11=u10.*loga(1:113,1) ;
u12=u10.^2-b(12,1)^2 ;
u13=xdata(2:114,9).^2-b(13,1)^2 ;
u14=xdata(2:114,7).^2-b(14,1)^2*xdata(2:114,9).^2 ;
u15=xdata(2:114,8).^2-b(15,1)^2*xdata(2:114,9).^2 ;
u16=xdata(2:114,10).^2-b(16,1)^2*xdata(2:114,9).^2 ;
u17=xdata(2:114,10).^2-b(17,1)^2*xdata(2:114,11).^2 ;

u=[ u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 ] ;
