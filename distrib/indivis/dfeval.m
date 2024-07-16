function df=dfeval(f,f0,x0,dx)

[n,nc]=size(f0) ;
[k,nc]=size(x0) ;

df=zeros(n,k) ;

if dx==0
   ax=abs(x0) ;
   dx=(max([ ax ones(k,1)*(1e-2) ]')')*(1e-5) ;
end

dxm=eye(k).*(dx*ones(1,k)) ;

for j=1:k
   dxj=dxm(:,j) ;
   df(:,j)=feval(f,x0+dxj)-feval(f,x0-dxj) ;
end

df=df./(ones(n,1)*(dx')*2) ;
