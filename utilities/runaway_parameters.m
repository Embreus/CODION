function [Ec, vc1, vc2] = runaway_parameters(params)

% Ec = Ec/E_D
% vc = vc/v_Ta

rhos = params.rhos; %sums to 1 <-> quasi-neutrality
Zs = params.Zs;
ms = params.ms; %in units of proton masses
Ts = params.Ts;

me = 9.10938291e-31 / 1.67262178e-27; 
ma = ms(1);
Za = Zs(1);
kappas = sqrt(ms*Ts(1)./(ma*Ts));
kappa_e = sqrt(me*Ts(1)/ma);

Zeff = sum(rhos.*Zs);
function f = F(x) 
    %calculates tau_ae/(m_a v_Ta) F_c
    G = @(r) (erf(r) - 2/sqrt(pi) * r.*exp(-r.^2))./(2*r.^2);
    f = 2*(1+me/ma) * G(kappa_e*x); %F = E^*/E_D
    for i=1:length(rhos)
        f = f + 2*Zs(i)*rhos(i)/Ts(i) * (1+ms(i)/ma) * G(kappas(i)*x);
    end 
end
% x=linspace(1e-5,100,10000);
% figure(10);clf;plot(x,F(x));

options = optimset('Display', 'off','LargeScale','off') ;
v_min = fminunc(@(x)abs(F(x)),15,options);
Ec = abs(Za/(2*Ts(1)*(1-Za/Zeff))*F(v_min));

vc1 = abs(fminunc(@(x)abs(F(x)-abs(params.EHat)),v_min-3,options));
vc2 = abs(fminunc(@(x)abs(F(x)-abs(params.EHat)),v_min+3,options));

end