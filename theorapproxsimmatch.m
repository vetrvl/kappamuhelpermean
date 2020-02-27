% kappa mu sim only

clear
close all
clc
Ps1=2; Ps2=1;
SNR=-35:2:2;
totalbits=1e5;
for kappa=[0,2,4]
    for mu=[1,3]
        jj=1;
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
        
        
        for N1=[8 16 32]
            N2=N1;
            for iii=1:length(SNR)
                sigma=sqrt(((Ps1+Ps2))/(4*(10^(SNR(iii)/10))));
                
                error=0;
                SNRlin=10.^(SNR(iii)/10);
                for iterations=1:totalbits
                    x1=2*randi([0,1])-1;
                    x2=2*randi([0,1])-1;
                    n=(randn + 1i*randn);
                    N0=(Ps1+Ps2)/(2*SNRlin);
                    A= sqrt(overallvariance)*randn+N1*overamean;  %         A=var1*randn+mean1;
                    B= sqrt(overallvariance)*randn+N2*overamean;
                    r =  A*sqrt(Ps1)*x1 +  B*sqrt(Ps2)*x2 + sigma*n;
                    
                    % Maximum Likelihood Detector
                    
                    % metrics=zeros(1,M);
                    
                    metrics11=norm(r - A*sqrt(Ps1) - B*sqrt(Ps2))^2;% 1 1 case
                    metrics00=norm(r + A*sqrt(Ps1) + B*sqrt(Ps2))^2; %-1 -1
                    metrics10=norm(r - A*sqrt(Ps1) + B*sqrt(Ps2))^2; % -1 1
                    metrics01=norm(r + A*sqrt(Ps1) - B*sqrt(Ps2))^2;% 1 -1
                    metrica=[metrics11,metrics00,metrics10,metrics01];
                    [val,indx]=min(metrica);
                    if indx ==1
                        decodebit=1;
                    elseif indx ==2
                        decodebit=1;
                    elseif indx ==3
                        decodebit=-1;
                    elseif indx ==4
                        decodebit=-1;
                    end
                    if x1*x2~=decodebit
                        error=error+1;
                    end
                end
                errorvec(iii)= error/totalbits;
                %         errorvecTheo(iii)=(qfunc(sqrt(Ps2)/sqrt(N0)));
                %         qfuncter(iii)=qfunc(sqrt(Ps2)/sqrt(N0)) - 0.5*qfunc((2*sqrt(Ps1)+sqrt(Ps2))/(sqrt(N0)))+0.5*qfunc((2*sqrt(Ps1)-sqrt(Ps2))/(sqrt(N0)));
                
                stringer=strcat("Simulation N1= N2= ",num2str(N1)," kappa= ", num2str(kappa),"mu= ",num2str(mu));
                
                sigma=sqrt(((Ps1+Ps2))/(4*(10^(SNR(iii)/10))));
                N0=sigma^2;
                NP=N0;
                %%% now add theory %%%%overamean
                a=(2*Ps2/(Ps1+Ps2));
                fun1 = @(t)  ((1./(1+ 2*N1*overallvariance.*a./ (4*(NP/Ps1)*((sin(t)).^2)) )).^(0.5)) .*...
                    exp(- ( (2*N1*overamean)^2 .* a ./(16*(NP/Ps1)*((sin(t)).^2))) ./...
                    (1+ 2*N1*overallvariance.*a./ (4*(NP/Ps1)*((sin(t)).^2)) ) ).*...
                    ((1./(1+ 2*N2*overallvariance.*a./ (4*(NP/Ps2)*((sin(t)).^2)) )).^(0.5)) .* ...
                    exp(- ( (2*N2*overamean)^2 .* a ./(16*(NP/Ps2)*((sin(t)).^2))) ./...
                    (1+ 2*N2*overallvariance.*a./ (4*(NP/Ps2)*((sin(t)).^2)) ) );
                peZ1(iii)=(1/pi)*integral(fun1,0,pi/2);
            end
            semilogy(SNR,peZ1,'-','linewidth',1.6,'DisplayName',"theory approx")
            
            hold on
            semilogy(SNR,errorvec,'o','linewidth',1.4,'DisplayName',stringer)
        end
    end
end
grid on
legend show
ylim([1e-6 1])
