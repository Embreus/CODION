function F = gen_distribution_function(f,theta,grid)

Ny = grid.Ny;
Nxi = grid.Nxi;
Nt = length(f(1,:));

Nth=length(theta);
F=zeros(Ny,Nth,Nt);

for th=1:Nth
F(1,th,:)=f(end,:);
end

P=zeros(Nxi,Nth);
for l=0:Nxi-1
    P0=legendre(l,cos(theta));
    P(l+1,:)=P0(1,:);
end

for tau=1:Nt
    %for Y=2:Ny
    for th = 1:Nth
        for L=0:Nxi-1
            F(2:Ny,th,tau)=F(2:Ny,th,tau)+f(1+L*(Ny-1):(Ny-1)+L*(Ny-1),tau)*P(L+1,th);
        end
    end
end
F = F/pi^(3/2);

end