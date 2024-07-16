'Christiano and Eichenbaum (1992) Stochastic Growth Case'
'Parameters: Theta,Alpha,Delta,lnGamm,sigx,aG,rhoG,sigG,6 Data Mom'

% SET INITIAL PARAMETER VALUES

bstart=[ 0.005 0.6 0.02 0.005 0.01 200 0.95 0.02 0.02 1.0 1.0 1.0 1.0 1.0]' ;

% SET GMM PARAMETERS

w0flag=0 ;

% w0flag=0 : W0=I is used as the initial weighing matrix
% w0flag=1 : bstart is used to calculate initial W0
% w0flag=2 : W0 passed to gmmest is used as initial W0.

w0=0 ; % initial W0, set to 0 if w0flag=0 or 1

nstep=2 ;

% Sets # of steps over weighting matrix, W0.
% Set w0flag=0 and nstep=2 to execute usual 2-Stage GMM.

miter=400 ;

% Sets maximum number of iterations on each GMM step

xtol=1e-6 ;

% gradient tolerance for convergence on each step

s0method=2;

% This variable is used to choose the method to calculate S0
% s0method=0  No lags used in computing s0
% s0method=1  Flat window used when computing s0 (may not be psd)
% s0method=2  Bartlett window used when computing s0
% s0method=3  Parzen window used when computing s0

lags=5 ;  % number of lags used in window methods for calculating s0

sstol=1e-12 ;

% step length tolerance - see quadmin

% LOAD THE DATA

load hoarding.dat ;

c=hoarding(:,1) ;
dk=hoarding(:,2) ;
g=hoarding(:,3) ;
y=hoarding(:,4) ;
n=hoarding(:,5) ;
k=hoarding(:,6) ;
hc=hpfilter(log(c),1600) ;
hi=hpfilter(log(dk),1600) ;
hg=hpfilter(log(g),1600) ;
hy=hpfilter(log(y),1600) ;
hn=hpfilter(log(n),1600) ;
hapl=hy-hn ;

global xdata ;

xdata=[ c dk g y n k hc hi hg hy hn hapl ] ;

% ESTIMATE THE MODEL, PERFORM TESTS, GRAPH IMPULSE RESPONSES

[b,v,q]=gmmest(bstart,nstep,miter,xtol,w0flag,w0,s0method,lags,sstol) ;

dmmy=cetests(b,v) ;

