% is gaussan square a nc-chisq matched except gaussian
clear
close all
mu=2;
va=4;
for mc=1:10000
    h(mc)=mu+sqrt(va)*randn;
    mag_h(mc)=h(mc)^2;
    ncchigen(mc)=va*random('Noncentral Chi-square',1,(mu^2)./(va));
    %% now cluster 2
    g(mc)=mu+sqrt(va)*randn;
    mag_g(mc)=g(mc)^2;
    overallmag(mc)= h(mc)^2 + g(mc)^2;
    overallchi(mc)=va*random('Noncentral Chi-square',2*1,2*(mu^2)./(va));
    rootoverallmag(mc)= sqrt(overallmag(mc)/2);
    ricipdf(mc)=sqrt(va/2)*random('Rician',sqrt(2*(mu^2)./(va)),1);
end
histogram(mag_h)
hold on
histogram(ncchigen)
figure
histogram(overallmag)
hold on
histogram(overallchi)
figure
histogram(rootoverallmag)
hold on
histogram(ricipdf)
