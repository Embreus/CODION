function [Ec, vc1, vc2] = runaway_parameters(params,settings)
rhos = params.rhos; %sums to 1 <-> quasi-neutrality
Zs = params.Zs;
ms = params.ms; %in units of proton masses
Ts = params.Ts;
me = 9.10938291e-31 / 1.67262178e-27; 
ma = ms(:,1);
Za = Zs(:,1);
EHat = abs(params.EHat);

switch settings.units
    case 'SI'
    eps0 = 8.85418782e-12; %m^-3 kg^-1 s^4 A^2, permittivity vacuum
    e_c  = 1.60217657e-19; %C, electron charge
    n_e  = -params.rhos(end);
    Zeff = Zs(1:end-1)*rhos(1:end-1)'/n_e;
    ln_Lambda = log(4*pi/3 * eps0^(3/2)/e_c^3 * (e_c*Ts(1))^(3/2)/n_e^(1/2));
    E_D  = ln_Lambda * n_e/(4*pi) * e_c^3/eps0^2 * 1/(e_c*Ts(end)); %V/m
    EHat =  2/Za * Ts(1)/Ts(end)*abs(1-Za/Zeff)*EHat/E_D;
    rhos = rhos/n_e;
    Ts   = Ts / Ts(end);
    otherwise
end  


kappas = sqrt(ms*Ts(1)./(ma*Ts));

if settings.electronCollisions
    Zeff = sum(rhos.*Zs);
    nbar = ms(1)*sum(rhos.*Zs./ms);
else
    Zeff = rhos(1:end-1)*Zs(1:end-1)';
    nbar = ms(1)*(rhos(1:end-1)./ms(1:end-1))*Zs(1:end-1)';
end
function f = F(x) 
    %calculates normalized collisional friction tau_ae/(m_a v_Ta) F_c, 
    %corresponding to E being normalized to EHat
    f = 0;
    G = @(r) (erf(r) - 2/sqrt(pi) * r.*exp(-r.^2))./(2*r.^2);
    if settings.electronCollisions
        kappa_e = sqrt(me*Ts(1)/ma);
        f = 2*(1+me/ma) * G(kappa_e*x); 
    end
    for i=1:length(rhos)
        f = f + 2*Zs(i)*rhos(i) * (1+ms(i)/ma) * G(kappas(i)*x);
    end 
end

options = optimset('Display', 'off','LargeScale','off') ;
v_min_g = (Zeff+nbar)^(1/3) * (9*pi/4 * ma/me)^(1/6);
v_min   = fminunc(@(x)abs(F(x)),v_min_g,options);
Ec = abs(Za/(2*Ts(1)*(1-Za/Zeff))*F(v_min)); %= Ec/E_D

vc1 = zeros(size(EHat));
vc2 = zeros(size(EHat));
for K = 1:length(EHat)
    if EHat(K) > F(v_min)
        %use solution to approximative Eq. as initial guess
        A = 4/(3*sqrt(pi)) * sqrt(Ts(1)*me/ma);
        x=sym('x');
        S = sort(real(double(solve(EHat*x^2-A*x^3-(nbar+Zeff) == 0))));
        vc1_g = S(2);
        vc2_g = S(3);
        %vc1_g = sqrt((Zeff+nbar)./EHat);
        %vc2_g = sqrt(ma/(me*Ts(1))) * 3*sqrt(pi)/4 .* EHat;
    else
        vc1_g = v_min;
        vc2_g = v_min;
    end    
    vc1(K) = abs(fminunc(@(x)abs(F(x)-EHat(K)),vc1_g(K),options));
    vc2(K) = abs(fminunc(@(x)abs(F(x)-EHat(K)),vc2_g(K),options));
end

end