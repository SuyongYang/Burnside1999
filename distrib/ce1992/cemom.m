function mm=cemom(b)

beta_=1.03^(-.25) ;
rho=[ 0  0 ; 0 b(7,1) ] ;
sigma=[ b(5,1) 0 ; 0 b(8,1) ] ;
sigma=sigma*sigma ;

drules=cesolve(beta_,b(1,1),b(2,1),b(4,1),b(3,1),b(6,1),rho) ;

indgro=[ -1 ; 1 ; 0 ; 1 ; 1 ; 1 ; 0 ; 0 ] ;
[facv,facr]=dshpmom(sigma,drules,indgro,0) ;

sd=sqrt(diag(facv)) ;

mm=[ sd(5,1)/sd(7,1) ; sd(9,1)/sd(7,1) ; sd(6,1)/sd(7,1) ;
     sd(6,1)/sd(8,1) ; sd(3,1)/sd(7,1) ; sd(7,1) ; facr(6,8) ] ;

