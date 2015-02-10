function [Ec, vc1, vc2] = runaway_parameters(params,settings)
rhos = params.rhos; %sums to 1 <-> quasi-neutrality
Zs = params.Zs;
ms = params.ms; %in units of proton masses
Ts = params.Ts;
if settings.units == 1
    rhos = rhos/abs(rhos(end));
    Ts = Ts / Ts(end);
end  
me = 9.10938291e-31 / 1.67262178e-27; 
ma = ms(1);
Za = Zs(1);
kappas = sqrt(ms*Ts(1)./(ma*Ts));

if settings.electronCollisions
    Zeff = sum(rhos.*Zs);
    nbar = ms(1)*sum(rhos.*Zs./ms);
else
    Zeff = rhos(1:end-1)*Zs(1:end-1)';
    nbar = ms(1)*(rhos(1:end-1)./ms(1:end-1))*Zs(1:end-1)';
end
function f = F(x) 
    %calculates tau_ae/(m_a v_Ta) F_c
    f=0;
    G = @(r) (erf(r) - 2/sqrt(pi) * r.*exp(-r.^2))./(2*r.^2);
    if settings.electronCollisions
        kappa_e = sqrt(me*Ts(1)/ma);
        f = 2*(1+me/ma) * G(kappa_e*x); %F = E^*/E_D
    end
    for i=1:length(rhos)
        f = f + 2*Zs(i)*rhos(i) * (1+ms(i)/ma) * G(kappas(i)*x);
    end 
end

options = optimset('Display', 'off','LargeScale','off') ;
v_min = fminunc(@(x)abs(F(x)),15,options);
Ec = abs(Za/(2*Ts(1)*(1-Za/Zeff))*F(v_min));

EHat = abs(params.EHat);
%use analytic approximative forms as initial guess
vc1_g = sqrt((Zeff+nbar)/EHat);
vc2_g = sqrt(ma/me)*(1/Ts(1))^(3/2) * 3*sqrt(pi)/4 * EHat;

vc1 = abs(fminunc(@(x)abs(F(x)-EHat),vc1_g,options));
vc2 = abs(fminunc(@(x)abs(F(x)-EHat),vc2_g,options));

end