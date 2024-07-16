function s0=s0calc(u,wt,wl)

% wt indicates the window type
%   0: no lags used in computing s0
%   1: flat window used-may not be psd
%   2: Bartlett window used
%   3: Parzen window used
% wl indicates the lag length for the window (largest non-zero lag)

[T,k]=size(u) ;
s0=(u'*u)/T ;
if wt~=0
   if wt==1
      for j=1:wl
         s0j=(u(1+j:T,:)'*u(1:T-j,:))/T ;
         s0=s0+s0j+(s0j') ;
      end
   elseif wt==2
      for j=1:wl
         s0j=(u(1+j:T,:)'*u(1:T-j,:))/T ;
         s0=s0+(s0j+(s0j'))*(1-j/(wl+1)) ;
      end
   else
      for j=1:wl
         s0j=(u(1+j:T,:)'*u(1:T-j,:))/T ;
         if j<((wl+1)/2)
            fac=1-6*(j/(wl+1))^2+6*(j/(wl+1))^3 ;
         else
            fac=2*((1-j/(wl+1))^3) ;
         end
         s0=s0+(s0j+(s0j'))*fac ;
      end
   end
end
