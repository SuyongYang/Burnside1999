function [x1,f1,niter]=quadmin(func,x0,xtol,nitermax,sstol) ;

global wmatrix_ ;

[k,nc]=size(x0) ;

g0=feval(func,x0) ;
dg0=dfeval(func,g0,x0,0) ;

f0=g0'*wmatrix_*g0 ;
df0=g0'*wmatrix_*dg0 ;

h0=inv(dg0'*wmatrix_*dg0) ;

xi=-h0*(df0') ;
convcrit=0 ;
niter=1 ;

while convcrit==0

  s=1 ;
  while s>sstol
    dx=s*xi ;
    x1=x0+dx ;
    g1=feval(func,x1) ;
    f1=g1'*wmatrix_*g1 ;
    if f1<f0
      g0=g1 ;
      f0=f1 ;
      x0=x1 ;
      s=0 ;
    end
    s=s*0.5 ;
  end
  test_=max( abs(xi)./( max([ abs(x0) ones(k,1) ]')' ) ) ;
  if niter==nitermax
    convcrit=1 ;
  elseif test_<xtol
    convcrit=1 ;
  else
    dg0=dfeval(func,g0,x0,0) ;
    df0=g0'*wmatrix_*dg0 ;
    h0=inv(dg0'*wmatrix_*dg0) ;
    xi=-h0*(df0') ;
  end
  niter=niter+1 ;
end
