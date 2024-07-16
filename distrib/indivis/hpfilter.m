function y=hpfilter(x,lamb)
[d,k]=size(x) ;
a=zeros(d,d) ;
for i=3:d-2
 a(i,i)=6*lamb+1;
 a(i,i+1)=-4*lamb;
 a(i,i+2)=lamb;
 a(i,i-2)=lamb;
 a(i,i-1)=-4*lamb;
end
a(2,2)=1+5*lamb;
a(2,3)=-4*lamb;
a(2,4)=lamb;
a(2,1)=-2*lamb;
a(1,1)=1+lamb;
a(1,2)=-2*lamb;
a(1,3)=lamb ;
a(d-1,d-1)=5*lamb+1;
a(d-1,d)=-2*lamb;
a(d-1,d-2)=-4*lamb;
a(d-1,d-3)=lamb;
a(d,d)=1+lamb;
a(d,d-1)=-2*lamb;
a(d,d-2)=lamb;
y=(eye(d)-inv(a))*x ;
