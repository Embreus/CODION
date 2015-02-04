%%% This script demonstrates how to specify parameters and settings and run
%%% CODION, obtaining the time-evolution of the ion distribution function 


rho_c = .3; %fraction of electrons due to carbon impurities
params.rhos = [1-rho_c rho_c]; %rhos_i = n_i * Z_i
params.Zs = [1 6]; %ion charges
params.ms = [2 12]; %ion masses, in units of proton masses. Only their ratios 
                    %appear. Need to add electrons with mass 1/1836 unless
                    %settings.electronCollisions = 1, which hard-codes it
params.Ts = [1 1]; %temperatures, Ts_i = T_i/T_e
params.nes = 1; %list of time-dependent electron densities. Don't touch

Zeff = params.Zs * params.rhos'; 

grid.Nxi = 30; %number of legendre modes used to represent the distribution
grid.Ny = 150; %number of velocity grid points
grid.yMax = 15; %maxium velocity on the grid
grid.Nt = 10; %number of time steps
grid.tMax = 1000; %final time, normalized to \tau_ie
tHat = linspace(0,grid.tMax,grid.Nt); %the list of times at which the
                                      %distribution will be evaluated


FHat = .08; %F/Za e E_D = E*/E_D = (1-Z/Zeff)E/E_D, effective electric field
EHat = 2/params.Zs(1) * params.Ts(1,1) * FHat; %normalized effective field
params.EHat = -EHat; %the normalized electric field EHat = 2/Z T/Te E*/E_D
                     %Silly sign convention of acceleration direction.
params.refresh_times = 0; % =/= 0 when time-dependent parameters

%%%% example of time-dependent electric field
%note that runaway_parameters.m does not work with time-dependent
%input, so gridMode = 'auto' can not be used
%params.refresh_times = linspace(0,grid.tMax,grid.Nt)';
%params.EHat = -EHat * sin( linspace(0,pi,grid.Nt) )';


settings.gridMode = 0; %can use 'auto' to automatically calculate 
                       %appropriate Ny and yMax to have well-converged
                       %solutions
settings.approximateCollOp = 0; %1 uses the v_Ti << v << v_Te asymptotic 
                                %form of the collision operator
settings.electronCollisions = 1; %automatically accounts for electron 
                                 %collisions; only need to specify ion 
                                 %composition in the params construct.
settings.momentumConservation = 1; %momentum-conserving self-collisions
settings.energyConservation = 1; %energy-conserving self-collisions
settings.energyConstant = 0; %1 removes energy from the system so that its 
                             %initial value is retained -- experimental
                             %setting, do not use.
settings.initialDistribution = 0; %Maxwellian initial distribution. 1 uses
                                  %a shifted Maxwellian.
settings.timeAdvanceMethod = 0; %0 is first-order backward differentiation 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,f] = CODION(grid,params,settings); %velocity grid x = v/v_Ta and 
                                      %distribution function f=f/(v_Ta^3 n_a)
grid.Ny = length(x);
grid.yMax = x(end); %just recalculates grid parameters in case of
                    %settings.gridMode = 'auto'

n = N_x1x2(x,f,0,1e4); %calculates number of particles with velocity 
                       %between 0 and 1e4 (i.e. total density)

Nth=200;
theta=linspace(0,pi,Nth); %creates a grid in theta on which the distribution
                          %is evaluated

fprintf('-----------------------------------\n')
fprintf('Generating distribution function...\n')
time=tic;
F = gen_distribution_function(f,theta,grid);
fprintf('Finished in %g seconds.\n',toc(time))

fprintf('Density conserved to %2.3g %%. \n',100*(n(end)-n(1))/n(1))


[Ec,vc1,vc2] = runaway_parameters(params); %numerically calculates Ec, vc1
                       %and vc2 from the full single-particle friction
                       %force.
n_RI = N_x1x2(x,f,vc1,1e4); %runaway density, density of particles with 
                            %v > vc1

AX = [-5 grid.yMax 1e-5  1];
figure(4)
set(gcf, 'Position', [10, 100, 600, 500],'color','w')
for tau=1:grid.Nt
    clf
    semilogy(x,abs(F(:,1,tau)),'k','linewidth',3)
    hold on
    semilogy(-x,abs(F(:,end,tau)),'k','linewidth',3)
    axis(AX)
    xlabel('v_{//} / v_{Ti}','fontsize',20,'fontweight','bold')
    ylabel('f_i / n_i','fontsize',20,'fontweight','bold')
    set(gca,'fontsize',16,'fontweight','bold','linewidth',3)
    pause(5e-2)
end


[K, TH] = meshgrid(x,theta);
figure(2)
clf
set(gcf, 'Position', [700, 100, 700, 400],'color','w')
r=linspace(-5,-1,20);
AX = [-5 17 -10 10];
lgF = log10(abs(F));
for tau=1:grid.Nt%1:1:grid.Nt
    clf
    contourf(K.*cos(TH),K.*sin(TH),lgF(:,:,tau)',r)
    hold on
    contourf(K.*cos(TH),-K.*sin(TH),lgF(:,:,tau)',r)
    colormap(flipud(GeriMap(100)))
    colorbar('ytick',-8:0,'YTickLabel',{'1e-8','1e-7','1e-6','1e-5','1e-4','1e-3', '1e-2', '1e-1','1'},'fontsize',16,'fontweight','bold','linewidth',2)
    caxis([r(1) r(end)])
    axis equal
    axis(AX)
    xlabel('v_{//} / v_{Ti}','fontsize',20,'fontweight','bold')
    ylabel('v_\perp / v_{Ti}','fontsize',20,'fontweight','bold')
    set(gca,'fontsize',16,'ytick',[-10 -5 0 5 10],'fontweight','bold','linewidth',2)
    pause(1e-4)
end




