function N = N_x1x2(x,f,x1,x2)
Ny = length(x);
Nt = length(f(1,:));

x1_ind = 1;
x2_ind = Ny;
for k = 1:Ny
    if x(k) > x1
        x1_ind = k-1;
        break
    end
end
for k = 1:Ny
    if x(k) > x2
        x2_ind = k;
        break
    end
end
ys = x(x1_ind:x2_ind);
ws = simpson_quad(ys);

F0 = zeros(Ny,Nt);
F0(1,:) = 4/sqrt(pi) * f(end,:);
F0(2:Ny,:) = 4/sqrt(pi) * f(1:Ny-1,:);
F = F0(x1_ind:x2_ind,:);

N = ws .* ys.^2 * F;
N = N';
end