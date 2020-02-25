% finding_meaner and variance using pdf
clc
clear all
close all
kappa=5;
mu=3; omega=1;
% alpha=0:0.01:3;
syms alpha
a=(1+kappa)^((mu+1)/2); 
b=alpha.^mu; 
c=mu*(1+kappa)*(alpha.^2)/omega; 
d=exp(-c); 
e=kappa^((mu-1)/2);
f=omega^((mu+1)/2);
g=2*mu*sqrt(kappa*(1+kappa)/omega); 
p=besseli(mu-1,g*alpha); 
pdf_e=(2.*mu.*a.*b.*d.*p)./(e.*f.*exp(mu*kappa));
% hold on; 
% plot(alpha, pdf_e);

meaner_theory=vpaintegral(alpha.*pdf_e,alpha,0,inf)
varia_theory=vpaintegral(alpha.*alpha.*pdf_e,alpha,0,inf)-meaner_theory.^2

for mc=1:1e5
    H(mc)=kappa_mu_channel(kappa,mu,1,1);
end
meaner=mean(abs(H))
varia=var(abs(H))
