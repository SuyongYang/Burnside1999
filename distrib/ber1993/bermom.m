function mm=bermom(b)

sigma=[ b(12,1) 0 ; 0 b(9,1) ] ; sigma=sigma*sigma ;
sg=exp(b(6,1)-b(4,1)) ;
rho=[ b(11,1) 0 ; 0 b(8,1) ] ;

drules=bersolve(b(1,1),b(2,1),b(3,1),b(5,1),sg,rho) ;

[facv,facr]=tshpmom(sigma,drules,0) ;

sd=sqrt(diag(facv)) ;

mm=[ sd(8,1) ; sd(6,1)/sd(8,1) ; sd(10,1)/sd(8,1) ; sd(2,1)/sd(8,1) ; sd(2,1)/sd(9,1) ] ;

