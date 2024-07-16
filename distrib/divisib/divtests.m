function tsts=divtests(b,varb)

mm=divmom(b) ;
[nm,nm2]=size(mm) ;
gm=dfeval('divmom',mm,b,0) ;

dm=divdata(b) ;
gd=dfeval('divdata',dm,b,0) ;

vm=zeros(nm,1) ;
vd=zeros(nm,1) ;

test=zeros(nm,1) ;
pv=zeros(nm,1) ;

for k=1:nm
  vm(k,1)=gm(k,:)*varb*(gm(k,:)') ;
  vd(k,1)=gd(k,:)*varb*(gd(k,:)') ;
  testn=mm(k,1)-dm(k,1) ;
  testd=(gm(k,:)-gd(k,:))*varb*((gm(k,:)-gd(k,:))') ;
  test(k,1)=(testn^2)/testd ;
  pv(k,1)=cdfchic(test(k,1),1) ;
end

sqvm=sqrt(vm) ;
sqvd=sqrt(vd) ;

'Model Moments   s.e.      Data Moments     s.e       Test      P-Value'

tsts=[mm sqvm dm sqvd test pv] ;
tsts 
