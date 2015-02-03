function ws = simpson_quad(x)
%array of quadrature weights ws,setting step length for last grid point
%to the value of the next last one
%int dx F  ->  ws' * F  
Ny = length(x);
dx = zeros(size(x));
dx(1:Ny-1) = x(2:end)-x(1:end-1);
dx(Ny) = dx(Ny-1);

ws = 48 * ones(size(x));
ws(1) = 17; ws(2) = 59; ws(3) = 43; ws(4) = 49;
ws((end-3):end) = ws(linspace(4,1,4));

ws =  ws .* dx / 48;

end