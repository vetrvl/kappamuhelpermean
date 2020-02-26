%%%%%

% gradstein table of integrals
% muu=2; ve=1; beta=1; alpha=3; % testing mu ve is now the problem





aa=gamma(muu+ve+0.5)./gamma(2.*ve);
bb= beta^(-1).*exp((beta.^(2))./(2.*alpha)).*alpha.^(-muu);
ww=whittakerM(-muu,ve,((beta.^(2))./(alpha)));
funcy=aa.*bb.*ww
% vpaintegral(funcy,,0,inf)





%% original funcy
syms x
xterm=x.^(muu-0.5);
expterm=exp(-alpha.*x);
Iterm=besseli(2.*ve,2.*beta.*sqrt(x));
funcyinteg=xterm.*expterm.*Iterm;

integ=2.*vpaintegral(funcyinteg,x,0,inf)
