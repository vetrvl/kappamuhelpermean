%% jth moment of kappa-mu SNR Matched
clear
clc
% mu=1;
kappa=4; jj=1;
mu=5;
% mu=(1+2*kappa)./(1+kappa.^2)
numerat=gamma(mu+(jj/2)).*exp(-kappa.*mu);
denom=gamma(mu).*((1+kappa).*mu).^(jj/2);
conflu=hypergeom(mu+(jj./2),mu,kappa.*mu);
overamean=(numerat./denom).*conflu;

jj=2;
numerat1=gamma(mu+(jj/2)).*exp(-kappa.*mu);
denom1=gamma(mu).*((1+kappa).*mu).^(jj/2);
conflu1=hypergeom(mu+(jj./2),mu,kappa.*mu);
overavar=(numerat1./denom1).*conflu1;
overallvariance=overavar-(overamean).^2;
%%%%%%%% now simula
N=10^4; 
m = sqrt( kappa/((kappa+1))) ; % Mean of in phase and quadrature 14 %component 15 
s = sqrt( 1/(2*(kappa+1)) ); %variance of in phase and quadrature 16 %component 17 
ni=zeros(1,N); nq=zeros(1,N);
for j=1:2*mu 
    if mod(j,2)==1
        norm1=m/sqrt(2)+s*randn(1,N); 
        ni=ni+norm1.^2;
    else
       norm2=m/sqrt(2)+s*randn(1,N);
       nq=nq+norm2.^2;
    end
end
h=sqrt(ni+nq)/sqrt(mu);
theta = 2*pi*rand(1,N);
h=h.*cos(theta)+sqrt(-1)*h.*sin(theta);
mean(abs(h))
overamean 
overallvariance
var(abs(h))
