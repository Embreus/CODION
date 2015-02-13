function F = gen_distribution_function(f,theta,grid)
Nx  = grid.Nx;
NL = grid.NL;
Nt  = length(f(1,:));
Nth = length(theta);
F   = zeros(Nx,Nth,Nt);

%Since f_l(0) = 0 for all l>0, CODION only stores the x = 0 value of
%f_0, at the end of the vector. Therefore, length(f(:,1)) = 1+(Ny-1)*Nxi
for th = 1:Nth
F(1,th,:) = f(end,:);
end

P = zeros(NL,Nth);
for l = 0:NL-1
    P0 = legendre(l,cos(theta));
    P(l+1,:) = P0(1,:);
end

for th = 1:Nth
    for L = 0:(NL-1)
        lowerLim = 1+L*(Nx-1);
        upperLim = (Nx-1)+L*(Nx-1);
        indicesForFL = lowerLim:upperLim;
        F(2:Nx,th,:) = squeeze(F(2:Nx,th,:)) + f(indicesForFL,:)*P(L+1,th);
    end
end
F = F/pi^(3/2); %raw CODION output normalized to F(0)=1

end