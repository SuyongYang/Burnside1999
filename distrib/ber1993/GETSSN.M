function ssn=getssn(theta_,alpha_,deltak,gammax,sg)

global beta_ capt xi ssf 

gammax=exp(gammax) ;
kyratio=(beta_*(1-alpha_))/(gammax-beta_*(1-deltak)) ;
si=(gammax-(1-deltak))*kyratio ;
sc=1-si-sg ;
ssn=alpha_/(theta_*sc*log(capt/(capt-xi-ssf))) ;
