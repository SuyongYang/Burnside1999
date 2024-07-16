function mm=indmom(b)

beta_=1.03^(-.25) ;
sigma=b(4,1) ; sigma=sigma*sigma ;

drules=indsolve(beta_,b(8,1),b(6,1),b(7,1),b(3,1)) ;

[facv,facr]=tshpmom(sigma,drules,0) ;

sd=sqrt(diag(facv)) ;

mm=[ sd(6,1) ; sd(4,1)/sd(6,1) ; sd(8,1)/sd(6,1) ; sd(5,1)/sd(6,1) ; sd(5,1)/sd(7,1) ; facr(5,7) ] ;
