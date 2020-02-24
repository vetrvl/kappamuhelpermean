% non zero mean abs complex gaussian
clc
clear 
totpower=1;
meane=0;
varia=1;
k=40;
s=sqrt(k/(k+1))*totpower;
sigma=totpower/sqrt(2*(k+1));
kappa1=s;mu1=1;
for mc=1:10000
% h(mc)=meane+sqrt(varia)*randn + 1j*randn;
% h(mc)=k+(randn + 1j*randn)./sqrt(2);
h(mc)=((sigma*randn+s)+1i*(randn*sigma+0)); %Rician Fading - single tap
magofnorm(mc)=abs(h(mc));
 Hhh(mc)=kappa_mu_channel(kappa1,mu1,1,1);
end
sigmasq=sigma.^2;
vsq=s.^2;
sigma.*sqrt(pi/2).*laguerreL(0.5,-vsq/(2*sigmasq))
mean(magofnorm)
% nowmixnagakami
mfac=(k.^2+2*k+1)/(2*k+1);
(gamma(mfac+0.5)./gamma(mfac)).*sqrt(totpower/mfac) 
