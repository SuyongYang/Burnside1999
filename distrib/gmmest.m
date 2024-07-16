function [b,v,q]=gmmest(b0,nstep,nmax,xtol,w0flag,w0,wtype,wlgth,sstol)

global wmatrix_ ;

[nb,nb1]=size(b0) ;
u0=gmmerr(b0) ;
[T,k]=size(u0) ;

if w0flag==0
   wmatrix_=eye(k) ;
elseif w0flag==1
   wmatrix_=inv(s0calc(u,wtype,wlgth)) ;
else
   wmatrix_=w0 ;
end

% 1st Step of GMM

[b,q,niter]=quadmin('mgmmerr',b0,xtol,nmax,sstol) ;

if niter>nmax
   'warning: maximum iterations reached'
end

u=gmmerr(b) ;
g0=mean(u)' ;
dg=dfeval('mgmmerr',g0,b,0) ;
idwd=inv(dg'*wmatrix_*dg) ;
s0=s0calc(u,wtype,wlgth) ;

'Iteration & Number of Steps'
[ 1 niter]

'Parameters'
[ b ]

iter=1 ;

while iter<nstep

 wmatrix_=inv(s0) ;
 b0=b ;

 [b,q,niter]=quadmin('mgmmerr',b0,xtol,nmax,sstol) ;

 if niter>nmax
  'warning: maximum iterations reached'
 end

 g0=feval('mgmmerr',b) ;
 dg=dfeval('mgmmerr',g0,b,0) ;
 v=(inv(dg'*wmatrix_*dg))/T ;
 sd=sqrt(diag(v)) ;

 iter=iter+1 ;

 'Iteration & Number of Steps'
 [iter niter]
 'GMM T*Q'
 q=T*q ;
 if k>nb
    [ q cdfchic(q,k-nb) ]
 else
    q
 end
 'Parameters & Std. Errors'
 [ b sd ]

end
