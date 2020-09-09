clear all ; close all ;clc ; 

gls= 0.00000; 
%work = (7.6*10^-6);
eps=1.0*10^-12 ;

max_eps=1.0*10^14 ; 

gls_kmin=7.6*10^-6 ; 
gls_cmu0=0.5477; 

gls_m=1.5 ;
gls_n=-1.0 ;
gls_p=3.0 ; 

ck=0.09 ;

gls_c2=1.92 ; 

cff3=1/gls_n; 
cff2=1.5 + (gls_m/gls_n);
cff1=3.0 + (gls_p/gls_n); 

tke=0.00000001; 

work=max(tke, gls_kmin); 

max_eps=1.0e14 ; 

dissip=(gls_cmu0.^cff1)*(work.^cff2)*(gls.^cff3) ;

%dissip= min(dissip, max_eps); 

inverse_dissip= max((1.0/dissip), max_eps) ; 

taufree=work*inverse_dissip; 

L=0.2; 

cff4=min( 1.0/( (ck.^2*tke).^(1.0/3.0)), max_eps); 

% CHANGED HERE

tauveg= ((L.^2)/(cff4));

taueff=min(taufree, tauveg) ; 

gls_veg=(gls_c2*tke/taueff)   

