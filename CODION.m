function OUT = CODION(grid0,params0,settings)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           CODION: COllisional Distribution of IONs
%           -------------------------------------------
%           Developed by Ola Embréus, 2014.
%           CODION paper: to be submitted to Physics of Plasmas
%           CODION MSc thesis: Ola Embréus 2014, electronically available 
%           at http://publications.lib.chalmers.se/records/...
%                                           fulltext/210276/210276.pdf
%
%           Discretization scheme originally written by 
%           Matt Landreman for CODE, September 2012.
%           CODE paper: M. Landreman et al. Comp. Phys. Comm. 185, 3 (2014)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('./utilities')

totalTime = tic;
switch settings.units 
    case 'SI'
        %Convert SI units to normalized values
        [grid,params,OUT] = NormalizeSIParameters(grid0,params0);
    otherwise
        grid   = grid0;
        params = params0;
end
settings.units = 0;
OUT.grid       = grid;
OUT.params     = params;
OUT.settings   = settings;

%this handles time-variable input, allowing only a
%few of the input parameters to vary, automatically 
%rescaling the other parameters to be of the same size.
[refresh_times,ms0,Zs0,Ts0,rhos0,nes0,EHat0] = RescaleParameters(params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% approximateCollOp    ; 1: Takes the v_Ti << v << v_Te limit of the CollOp 
%                           and neglects electron scattering, as is done in  
%                           the analytic calculation of f_RI.
% electronCollisions   ; 1: Automatically includes electron scattering.
% momentumConservation ; 1: Momentum conserving self-collision operator.
% energyConservation   ; 1: Energy conserving self-collision operator.
% initialDistribution  ; 0: Maxwellian, 1: Shifted Maxwellian (by x0)
% timeAdvanceMethod    ; 0: Implicit Euler
%                        1: BDF2 (Backward Differentiation Formula)
%                        2: Trapezoid rule (Adams-Moulton 2)
% gridMode             ; 0: Uniform grid with Ny, yMax and Nxi as given
%                   'auto': Uniform grid with Ny, yMax, Nxi automatically 
%                           calculated based on tMax to yield well-
%                           converged solution. Does not affect Nt.
% units              ;'SI': Input parameters are assumed to be in 
%                           dimensional units - T in eV, densities in m-3,
%                           EHat = electric field in V/m, tMax in seconds,
%                           mass in proton masses
%                otherwise: Assumes normalized input parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialDistribution
% 0 = A Maxwellian.
% 1 = A shifted Maxwellian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%timeAdvanceMethod
% 0 = Backward Euler = BDF1 = Adams-Moulton 1
% 1 = BDF2 (backward differentiation formula)
% 2 = Trapezoid rule ( = Adams-Moulton 2)

% dt0: timestep:
% tMax0: Duration of simulation:
% Nxi0 = number of Legendre polynomials in xi = v_|| / v:
% Ny0 = number of grid points in x = v / v_Ta:
% yMax0 = maximum value of y retained;

tMax  = grid.tMax;
Nt    = grid.Nt;
dt    = tMax/(Nt-1);
tHat  = linspace(0,tMax,Nt);
Nx    = grid.Nx;
xMax  = grid.xMax;
NL    = grid.NL;

% Boundary condition at yMax.
xMaxBoundaryCondition = 4;
% 1 = Dirichlet: F=0.
% 2 = Robin: dF/dy + (2/y)*F=0, which forces F to behave like 1/y^2.
% 3 = Do not apply a boundary condition at yMax.
% 4 = Dirichlet: F=0, with a bit of extra d2dy2 added at the last few grid
%       points to eliminate ringing.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of input parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('***************************************************************\n')
fprintf('Time-dependent CODION, using ')
switch settings.timeAdvanceMethod
    case 0
        % Backward Euler
        fprintf('backward Euler')
    case 1
        % BDF2
        fprintf('BDF2')
    case 2
        % Trapezoid rule
        fprintf('trapezoid rule')
    otherwise
        error('Invalid timeAdvanceMethod')
end
fprintf(' for time advance.\n')

%Print which settings for self-collision operator is used
fprintf('Using ')
if ~(settings.momentumConservation||settings.energyConservation)
    fprintf('only test-particle ')
end
if settings.momentumConservation
    fprintf('momentum conserving ')
end
if settings.energyConservation
    fprintf('and energy conserving ')
end
fprintf('self-collision operator.\n')

if settings.electronCollisions
    Zeff = params.Zs(1,:)*params.rhos(1,:)';
    nbar = params.ms(1)*(params.rhos./params.ms)*params.Zs';
else
    Zeff = params.Zs(1,1:end-1)*params.rhos(1,1:end-1)';
    nbar = params.ms(1)*(params.rhos(1,1:end-1)./params.ms(1,1:end-1))*params.Zs(1,1:end-1)';
end
OUT.Zeff = Zeff;
OUT.nbar = nbar;
params_tmp      = params;
params_tmp.EHat = max(abs(params.EHat));
[Ec,xc1,xc2]    = runaway_parameters(params_tmp,settings); %numerically calculates 
                      %Ec/E_D, vc1/v_Ta and vc2/v_Ta from the full single-particle 
                      %friction force.
OUT.Ec  = Ec;
OUT.vc1 = xc1;
OUT.vc2 = xc2;

EcHat   = Ec/abs(params.Zs(1)/(2*params.Ts(1)*(1-params.Zs(1)/Zeff)));
%print a couple of relevant physical parameters
fprintf('EHat: %2.3g, E/E_c: %2.3g, Zeff: %2.3g \n',EHat0(1), abs(EHat0(1)/EcHat), Zeff)
fprintf('xc1: %2.3g, xc2: %2.3g, nbar: %2.3g \n', xc1, xc2, nbar)


% Generate differentiation matrices.
xMin=0;
scheme = 12; %12: uniform grid
switch settings.gridMode
    case 0
        [x,~, ddx, d2dx2] = m20121125_04_DifferentiationMatricesForUniformGrid(Nx, xMin, xMax, scheme);
    case 1 % avoid using this setting, experimental non-uniform grid
        scheme = 12;
        [s,~, dds, d2ds2] = m20121125_04_DifferentiationMatricesForUniformGrid(Nx, xMin, 1, scheme);
        
        c0 = 1;
        c1 = 0.01;
        c2 = 5;
        x  = c1*s.^c2 + c0*s;
        dyds   = c2*c1*s.^(c2-1) + c0;
        d2yds2 = c2*(c2-1)*c1*s.^(c2-2);
        
        ddx    = diag(1./dyds)*dds;
        d2dx2  = -diag(d2yds2 ./ (dyds.^3)) * dds + diag((1./dyds).^2)*d2ds2;
        d2dx2(end,:) = 0;
    case 'auto' %overrides Nx, xMax and NL to set ''suitable'' values based on the values 
                %chosen for tMax and EHat, so that solutions will be well converged
        xMax = xc2+8;
        dx0  = 0.45;
        dx   = dx0/tMax^(1/4);
        Nx   = min([round(xMax/dx),1500]);
        NL0 = 40; %suitable option when xc2 \sim 10
        NL  = min([round(NL0*xc2/10),500]); %assuming constant width of runaway bump 
                                             %at accumulation point, Nxi \propto xc2
        [x,~, ddx, d2dx2] = m20121125_04_DifferentiationMatricesForUniformGrid(Nx, xMin, xMax, scheme);
        OUT.grid.NL = NL;
        OUT.grid.Nx = Nx;
        OUT.grid.xMax = xMax;
end
%print grid parameters
fprintf('Nx: %d,   xMax: %2.4g,   NL: %d,   dt: %2.3g,  tMax: %2.4g, Nt: %d\n',Nx,xMax,NL,dt,tMax,Nt)

% Make x a row vector:
x = x';
OUT.x = x;

% Order of rows in the matrix and right-hand side:
% --------------------
% for L=0:(NL-1)
%   for ix = 1:(Nx-2)
%     Impose kinetic equation
%   Impose boundary condition at xMax
% Enforce regularity: dF/dx=0 at x=0

% Order of columns in the matrix, corresponding to rows in the solution vector:
% --------------------
% for L=0:(NL-1)
%   for ix = 1:(Nx-1)
%     Value of F
% Value of F at y=0

matrixSize = NL*(Nx-1) + 1;

% Predict roughly how many nonzero elements will be in the sparse
% matrix. This speeds up the code by eliminating the need to reallocate
% memory during matrix construction.
% It is sensitive to which self-collision operator is used -- the 
% conserving ones add dense blocks.
predictedFillFactor = (3*NL*nnz(abs(ddx) + abs(d2dx2)))/(matrixSize*matrixSize);
if settings.momentumConservation && settings.energyConservation
    predictedFillFactor = predictedFillFactor + 2*Nx*Nx/(matrixSize*matrixSize);
elseif settings.momentumConservation || settings.energyConservation
    predictedFillFactor = predictedFillFactor + Nx*Nx/(matrixSize*matrixSize);
end

fprintf('Matrix size: %g\n',matrixSize)

OUT.f = zeros(matrixSize, Nt);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial condition:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fMinus1 = zeros(matrixSize,1);
switch settings.initialDistribution
    case 0
        % A stationary Maxwellian. This sets our normalization of f
        fMinus1(1:(Nx-2)) = exp(-x(2:(Nx-1)).^2);
        fMinus1(end) = 1;
    case 1
        % A shifted (non-normalized) Maxwellian
        x0 = .3; %flow velocity
        xVec = x(2:Nx-1);
        fMinus1(1:Nx-2)    = (1-x0^2).*exp(-xVec.^2);
        fMinus1(Nx:2*Nx-3) = 2 * x0 * xVec.*exp(-xVec.^2);
        fMinus1(end) = 1;
    case 'input'
        % Initializes with some distribution params.initialDistFunc = f;
        fMinus1 = params.initialDistFunc;
    otherwise
        error('Invalid setting for initial distribution')
end  
fMinus2 = fMinus1;
OUT.f(:,1)  = fMinus1;

ne0 = nes0(1,1);
Ta0 = Ts0(1,1);
ma0 = ms0(1,1);
me  = 1/1822.89;
kappa_e = sqrt(me*Ta0/ma0);

Phi = @(r) erf(r); %Error function = 2/sqrt(pi) int_0^r ds exp(-s^2)
dPhi = @(r) 2/sqrt(pi) * exp(-r.^2); 
G = @(r)(Phi(r)-r.*dPhi(r))./(2*r.^2); %Chandrasekhar function
dG =@(r) 2/sqrt(pi) * (1+1./(r.^2)).*exp(-r.^2) - Phi(r)./(r.^3);

x2 = x.*x;

refresh_counter = 0;
for iteration = 2:Nt
    if length(refresh_times) > refresh_counter  %check if it's time to 
                                                %refresh the matrix.
        if refresh_times(refresh_counter + 1) <= tHat(iteration)
            refresh_counter = refresh_counter + 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Begin building matrix.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ms   = ms0(refresh_counter,:);
            Ts   = Ts0(refresh_counter,:);
            rhos = rhos0(refresh_counter,:);
            Zs   = Zs0(refresh_counter,:);
            ne   = nes0(refresh_counter);
            EHat = EHat0(refresh_counter);
            
            kappas    = sqrt(ms*Ta0./(ma0*Ts));
            lnLambda0 = 1;
            lnLambda  = 1;
            collision_time_mod = ne * lnLambda / (ne0 * lnLambda0);

            % I'm sorry that this is a bit ugly -- lots of if statements in
            % this part of the code due to various settings, could possibly
            % be compactified
            if settings.electronCollisions
                coefficientOf_f = 2*Ta0*(2./x.*G(kappa_e*x) + kappa_e*dG(kappa_e*x));
                coefficientOf_dfdx = 2*Ta0*G(kappa_e*x) + G(kappa_e*x)./x2 + kappa_e*dG(kappa_e*x)./x;
                coefficientOf_d2fdx2 = G(kappa_e*x)./x;
            else
                coefficientOf_f = 0;
                coefficientOf_dfdx = 0;
                coefficientOf_d2fdx2 = 0;
            end
               
            
            if settings.approximateCollOp
                coefficientOf_f = zeros(size(x));
                coefficientOf_dfdx = (nbar./x.^2 - nbar./(2*x.^4));
                coefficientOf_d2fdx2 = nbar./(2*x.^3);
            else
                for i=1:length(rhos)
                    coefficientOf_f = coefficientOf_f + rhos(i)*Zs(i)*2*Ta0/Ts(i)*(2*G(kappas(i)*x)./x + kappas(i)*dG(kappas(i)*x));
                    coefficientOf_dfdx = coefficientOf_dfdx + rhos(i)*Zs(i)*(2*Ta0/Ts(i)*G(kappas(i)*x) + G(kappas(i)*x)./x2 + kappas(i)*dG(kappas(i)*x)./x);
                    coefficientOf_d2fdx2 = coefficientOf_d2fdx2 + rhos(i)*Zs(i)*G(kappas(i)*x)./x;
                end
            end
           

            energyScattering = (diag(coefficientOf_f) ...
                + diag(coefficientOf_d2fdx2)*d2dx2 + diag(coefficientOf_dfdx)*ddx);
            if xMaxBoundaryCondition==4
                fakeViscosity = exp((x-xMax)/(0.1));
                energyScattering = energyScattering + (1e-2)*diag(fakeViscosity)*d2dx2;
            end

            if settings.electronCollisions
                xPartOfPitchAngleScattering = (Phi(kappa_e*x)-G(kappa_e*x))./(2*x.^3);
            else
                xPartOfPitchAngleScattering = 0;
            end
            
            
            if settings.approximateCollOp
                xPartOfPitchAngleScattering = sum(rhos.*Zs)./(2*x.^3);
            else
                for i=1:length(rhos)
                   xPartOfPitchAngleScattering = xPartOfPitchAngleScattering + rhos(i)*Zs(i) * (Phi(kappas(i)*x) - G(kappas(i)*x))./(2*x.^3);
                end
            end
            
            xPartOfPitchAngleScatteringMatrix = diag(xPartOfPitchAngleScattering);

            energyScattering = collision_time_mod * energyScattering;
            xPartOfPitchAngleScatteringMatrix = collision_time_mod * xPartOfPitchAngleScatteringMatrix;

  
            %%%%% MOMENTUM AND ENERGY CONSERVATION (probably) NOT YET %%%%%
            %%%%% PROPERLY NORMALIZED WITH TIME-DEPENDENT PARAMETERS  %%%%%
            energyConservingLittleMatrix = zeros(Nx,Nx);
            if settings.energyConservation
                energyConservingLittleMatrix = GenerateEnergyConservingLittleMatrix(Ta0,Ts,Phi,dPhi,collision_time_mod,rhos,Zs,x);
            end
            momentumConservingLittleMatrix = zeros(Nx,Nx);
             if settings.momentumConservation
                momentumConservingLittleMatrix = GenerateMomentumConservingLittleMatrix(Ta0,Ts,Phi,dPhi,collision_time_mod,rhos,Zs,x);
            end


            % Initialize arrays for building the sparse matrix:
            sparseCreatorIndex = 1;
            estimated_nnz   = 0;
            sparseCreator_i = 0;
            sparseCreator_j = 0;
            sparseCreator_s = 0;
            resetSparseCreator(matrixSize,predictedFillFactor)

            if xMaxBoundaryCondition == 3
                rowRange = 1:(Nx-1);
            else
                rowRange = 1:(Nx-2);
            end

            for L=0:(NL-1)
                % Add collision operator
                xPartMatrix = energyScattering - L*(L+1) * xPartOfPitchAngleScatteringMatrix;
                if L==0
                    xPartMatrix = xPartMatrix + energyConservingLittleMatrix;
                end
                if L==1
                    xPartMatrix = xPartMatrix + momentumConservingLittleMatrix;
                end

                rowIndices    = L*(Nx-1) + rowRange;
                columnIndices = L*(Nx-1) + (1:(Nx-1));
                addSparseBlock(rowIndices, columnIndices, xPartMatrix(1+rowRange, 2:Nx))
                if L==0
                    addSparseBlock(rowIndices, matrixSize, xPartMatrix(1+rowRange, 1))
                end


                % Add electric field term- 
                % Sub-diagonal term in L:
                ell = L-1;
                if ell >=0
                    columnIndices = ell*(Nx-1) + (1:(Nx-1));
                    littleMatrix  = L/(2*L-1)* EHat*ddx  - diag(EHat*(L-1)*L/(2*L-1)./x);
                    addSparseBlock(rowIndices, columnIndices, littleMatrix(1+rowRange, 2:Nx))
                    if ell==0
                        addSparseBlock(rowIndices, matrixSize, littleMatrix(1+rowRange, 1))
                    end
                end

                % Add electric field term- 
                % Super-diagonal term in L:
                ell = L+1;
                if ell<NL
                    columnIndices = ell*(Nx-1) + (1:(Nx-1));
                    littleMatrix = (L+1)/(2*L+3)*EHat*ddx  + diag(EHat*(L+1)*(L+2)/(2*L+3)./x);
                    addSparseBlock(rowIndices, columnIndices, littleMatrix(1+rowRange, 2:Nx))
                end

            end

            operator = createSparse();
            resetSparseCreator(matrixSize,predictedFillFactor)

            indices = 1:(matrixSize-1);
            addToSparse(indices, indices, ones(size(indices)))

            % For the special point at x=0, apply Neumann condition dF/dx=0:
            addToSparse(matrixSize, matrixSize, ddx(1,1))
            addSparseBlock(matrixSize, 1:(Nx-1), ddx(1, 2:Nx))

            % Impose boundary conditions at x=xMax:
            if xMaxBoundaryCondition==2
                % Add Robin boundary condition dF/dx + (2/x)*F = 0 at xMax:
                for L=0:(NL-1)
                    rowIndex = L*(Nx-1) + Nx-1;
                    columnIndices = L*(Nx-1) + (1:(Nx-1));
                    boundaryCondition = ddx(Nx,:);
                    boundaryCondition(Nx) = boundaryCondition(Nx) + 2/xMax - 1; %Subtract 1 since we already put a 1 on the diagonal.
                    addSparseBlock(rowIndex, columnIndices, boundaryCondition(2:Nx))
                    if L==0
                        addSparseBlock(rowIndex, matrixSize, boundaryCondition(1))
                    end
                end
            end

            diagonalPlusBoundaryCondition = createSparse();
            
            switch settings.timeAdvanceMethod
                case 0
                    % Backward Euler
                    timeAdvanceMatrix = diagonalPlusBoundaryCondition - dt*operator;
                case 1
                    % BDF2
                    timeAdvanceMatrix = diagonalPlusBoundaryCondition - (2/3)*dt*operator;
                case 2
                    % Trapezoid rule
                    timeAdvanceMatrix = diagonalPlusBoundaryCondition - (dt/2)*operator;
                otherwise
                    error('Invalid timeAdvanceMethod')
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % End of building the matrix.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [factor_L, factor_U, factor_P, factor_Q] = lu(timeAdvanceMatrix);    
        end
    end %end of if statements determining whether matrix should be rebuilt
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Beginning of time-advance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Build the right-hand side
    switch settings.timeAdvanceMethod
        case 0
            rhs = fMinus1;
        case 1
            rhs = (4/3)*fMinus1 - (1/3)*fMinus2;
        case 2
            rhs = fMinus1 + (dt/2)*operator*fMinus1;
    end

    % Handle boundary condition at x=0:
    rhs(end) = 0;

    % Handle boundary condition at x = xMax:
    if xMaxBoundaryCondition ~= 3
        rhs((1:NL)*(Nx-1)) = 0;
    end

    % Now step forward in time.
    % The next line is equivalent to 'soln = matrix \ rhs', but much faster:
    soln = factor_Q * (factor_U \ (factor_L \ (factor_P * rhs)));
    OUT.f(:,iteration) = soln;

    fMinus2 = fMinus1;
    fMinus1 = soln;
    
end

fprintf('Ion distribution obtained in: %2.4gs \n',toc(totalTime))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Below are utility functions used in the main body of the code above %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: make the interface clear: which variables are written two
function resetSparseCreator(matrixSize,predictedFillFactor)
    sparseCreatorIndex = 1;
    estimated_nnz      = floor(matrixSize*matrixSize*predictedFillFactor);
    sparseCreator_i    = zeros(estimated_nnz,1);
    sparseCreator_j    = zeros(estimated_nnz,1);
    sparseCreator_s    = zeros(estimated_nnz,1);
end

function addToSparse(i,j,s)
    % Adds values to the sparse matrix.
    n=numel(i);
    if n ~= numel(j)
        error('Error A');
    end
    if n ~= numel(s)
        error('Error B');
    end
    if any(i<1)
        error('Error Q: i<1');
    end
    if any(j<1)
        error('Error Q: j<1');
    end
    sparseCreator_i(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = i;
    sparseCreator_j(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = j;
    sparseCreator_s(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = s;
    sparseCreatorIndex = sparseCreatorIndex+n;
    if sparseCreatorIndex > estimated_nnz
        fprintf('Warning! estimated_nnz is too small. Increase predictedFillFactor.\n')
    end
end

function addSparseBlock(rowIndices, colIndices, block)
    % Adds a block to the sparse matrix.
    % rowIndices and colIndices should be vectors.
    % numel(rowIndices) should equal the number of rows in 'block'.
    % numel(colIndices) should equal the number of columns in 'block'.
    s=size(block);
    if (s(1) ~= numel(rowIndices)) || (s(2) ~= numel(colIndices))
        error('Error in addSparseBlock!  Input sizes are not consistent.')
    end
    [rows, cols, values] = find(block);
    addToSparse(rowIndices(rows),colIndices(cols),values)
end

function sparseMatrix = createSparse()
    % After you are done adding elements to the sparse matrix using
    % addToSparse() and addSparseBlock(), call this function to
    % finalize the matrix.
    sparseMatrix = sparse(sparseCreator_i(1:(sparseCreatorIndex-1)), sparseCreator_j(1:(sparseCreatorIndex-1)), sparseCreator_s(1:(sparseCreatorIndex-1)), matrixSize, matrixSize);
    resetSparseCreator(matrixSize,predictedFillFactor)
end


function [refresh_times,ms0,Zs0,Ts0,rhos0,nes0,EHat0]=RescaleParameters(params)
    refresh_times = params.refresh_times; 
    ms0   = RescaleX(params.ms, refresh_times);
    Zs0   = RescaleX(params.Zs, refresh_times);
    Ts0   = RescaleX(params.Ts, refresh_times);
    rhos0 = RescaleX(params.rhos, refresh_times);
    nes0  = RescaleX(params.nes, refresh_times);
    EHat0 = RescaleX(params.EHat, refresh_times);
end

function A = RescaleX(X,refresh_times)
    N_rfrs = length(refresh_times);
    [N_rows, N_cols] = size(X);
    if N_rows < N_rfrs
        A_temp = zeros(N_rfrs,N_cols);
        A_temp(1:N_rows,:) = X;
        for k = N_rows+1:N_rfrs
            A_temp(k,:) = X(N_rows,:);
        end
        A = A_temp;
    else
        A = X;
    end
end


% TODO: package the input variables into one record to reduce the number of arguments
function MC_little_matrix = GenerateMomentumConservingLittleMatrix(Ta0, Ts, Phi, dPhi, collision_time_mod, rhos, Zs, x)
    fM     = exp( - x'.^2 * Ta0/Ts(1));
    ws     = simpson_quad(x);
    U_DOWN = ws .* x .* ( Phi(x) - x.*dPhi(x) ) * fM;
    Uterm_cols = ws .* ( Phi(x) - x.*dPhi(x) ) / U_DOWN;
    Uterm_rows = collision_time_mod * rhos(1) * Zs(1) * 4*G(x) .* fM';
    MC_little_matrix = Uterm_rows'*Uterm_cols;
end
    
% TODO: package the input variables into one record to reduce the number of arguments
function EC_little_matrix = GenerateEnergyConservingLittleMatrix(Ta0, Ts, Phi, dPhi, collision_time_mod, rhos, Zs, x)
    fM     = exp( - x'.^2 * Ta0/Ts(1));
    ws     = simpson_quad(x);
    Q_DOWN = ws .* (2*x.^3) .* ( Phi(x) - 2*x.*dPhi(x) ) * fM;
    Qterm_cols = ws .* (2*x) .* ( Phi(x) - 2*x.*dPhi(x) ) / Q_DOWN;
    Qterm_rows = collision_time_mod * rhos(1) * Zs(1) * 2./x .*(Phi(x) - 2*x.*dPhi(x)).* fM';
    Qterm_rows(1) = collision_time_mod * rhos(1) * Zs(1) * (-4/sqrt(pi));
    
    EC_little_matrix = Qterm_rows'*Qterm_cols;
end

function [grid,params,OUT] = NormalizeSIParameters(grid,params)
    %Changes parameters to CODION-normalization. Incompatibable 
    %with settings.electronCollisions = 1
    eps0 = 8.85418782e-12; %m^-3 kg^-1 s^4 A^2, permittivity vacuum
    m_p  = 1.67262178e-27; %kg, proton mass
    m_e  = 9.10938291e-31;
    e_c  = 1.60217657e-19; %C, electron charge
    
    n_e   = -params.rhos(end);
    m_a   = m_p*params.ms(1); %kg
    v_Te  = sqrt(2*e_c*params.Ts(end)/m_e);
    v_Ta  = sqrt(2*e_c*params.Ts(1)/m_a); %m/s
    ln_Lambda = log(4*pi/3 * eps0^(3/2)/e_c^3 * (e_c*params.Ts(1))^(3/2)/n_e^(1/2));
    nu_ie = ln_Lambda * n_e/(4*pi) * (params.Zs(1)*e_c^2/(m_a*eps0))^2 /v_Ta^3; %s^-1
    E_D   = ln_Lambda * n_e/(4*pi) * e_c^3/eps0^2 * 1/(e_c*params.Ts(end)); %V/m
    
    OUT.lnLambda = ln_Lambda;
    OUT.v_Te  = v_Te;
    OUT.v_Ta  = v_Ta;
    OUT.nu_ie = nu_ie;
    OUT.E_D   = E_D;
    
    grid.tMax   = grid.tMax * nu_ie;
    params.Ts   = params.Ts / params.Ts(end);
    params.rhos = params.rhos / abs(params.rhos(end));
    Zeff = params.Zs(1:end-1)*params.rhos(1:end-1)';
    params.EHat = (1-params.Zs(1)/Zeff)*params.EHat/E_D ...
        *2/params.Zs(1) * params.Ts(1)/params.Ts(end);
end

end