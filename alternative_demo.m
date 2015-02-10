%%% The script demonstrates another set of options,
%%% for example using SI units as input and automatic grid choice

clear
clc

nD = 2.42e19;   % m^-3
nC = 0.04*nD; 
ne = nD + nC*6;

me = 9.10938291e-31 / 1.67262178e-27; 
params.rhos = [nD nC*6 -ne];
params.Zs   = [1 6 -1];
params.ms   = [2 12 me];   %in units of proton masses
params.Ts   = [1 1 1]*1e3; %eV
params.nes  = 1;

Zeff = params.Zs(1:end-1) * params.rhos(1:end-1)'/ne;
nbar = params.ms(1)*sum((params.rhos(1:end-1).*params.Zs(1:end-1))./params.ms(1:end-1))/ne;

grid.Nxi  = 30;
grid.Ny   = [];
grid.yMax = [];
grid.Nt   = 20;
grid.tMax = 3; %seconds

tHat = linspace(0,grid.tMax,grid.Nt);

%params.refresh_times = linspace(0,grid.tMax,grid.Nt);
%params.EHat = EHat * sin(2*pi*params.refresh_times'/grid.tMax * 3) ;
params.refresh_times = 0;
params.EHat = -1.64; %V/m

settings.gridMode             = 'auto';
settings.approximateCollOp    = 0;
settings.electronCollisions   = 0;
settings.momentumConservation = 1;
settings.energyConservation   = 1;
settings.initialDistribution  = 0;
settings.timeAdvanceMethod    = 0;
settings.units                = 1; %1 = SI units

[x,f]     = CODION(grid,params,settings);
grid.Ny   = length(x);
grid.yMax = x(end);

[Ec,vc1,vc2] = runaway_parameters(params,settings);
n    = N_x1x2(x,f,0,1e4);
n_RI = N_x1x2(x,f,vc1,1e3);




Nth=200;
theta=linspace(0,pi,Nth);
fprintf('-----------------------------------\n')
fprintf('Generating distribution function...\n')
time=tic;
F = gen_distribution_function(f,theta,grid);
fprintf('Finished in %gs seconds.\n',toc(time))


pauseTime = 0.5;

AX = [-5 grid.yMax 1e-5  2e-1];
figure(1)
set(gcf, 'Position', [10, 100, 600, 500])
for tau=1:grid.Nt
    clf
    semilogy(x,abs(F(:,1,tau)),'k','linewidth',3)
    hold on
    semilogy(-x,abs(F(:,end,tau)),'k','linewidth',3)
    axis(AX)
    xlabel('v_{//} / v_{Ti}','fontsize',20,'fontweight','bold')
    ylabel('f_i / n_i','fontsize',20,'fontweight','bold')
    set(gca,'fontsize',16,'fontweight','bold','linewidth',3)
    pause(pauseTime)
end


[K, TH] = meshgrid(x,theta);
figure(4)
clf
set(gcf, 'Position', [700, 100, 600, 500])
r=linspace(-6,-3,15);
AX = [-5 15 -10 10];
lgF = log10(abs(F));
for tau=1:1:grid.Nt
    clf
    contourf(K.*cos(TH),K.*sin(TH),lgF(:,:,tau)',r)
    hold on
    colormap(flipud(GeriMap(100)))
    contourf(K.*cos(TH),-K.*sin(TH),lgF(:,:,tau)',r)
    caxis([r(1) r(end)])
    colorbar('ytick',-8:0,'YTickLabel',{'1e-8','1e-7','1e-6','1e-5','1e-4','1e-3', '1e-2', '1e-1','1'},'fontsize',16,'fontweight','bold','linewidth',2)
    axis equal
    axis(AX)
    xlabel('v_{//} / v_{Ti}','fontsize',20,'fontweight','bold')
    ylabel('v_\perp / v_{Ti}','fontsize',20,'fontweight','bold')
    set(gca,'fontsize',16,'ytick',[-10 -5 0 5 10],'fontweight','bold','linewidth',2)
    pause(pauseTime)
end

