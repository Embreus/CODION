function [x,f] = CODION(grid,params,settings)
addpath('./utilities')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           CODION: COllisional Distribution of IONs
%           -------------------------------------------
%           Developed by Ola Embréus, 2014.
%           CODION paper: to be submitted to Physics of Plasmas
%           CODION MSc thesis: Ola Embréus 2014, electronically available 
% at http://publications.lib.chalmers.se/records/fulltext/210276/210276.pdf
%
%           Discretization scheme originally
%           written by Matt Landreman for CODE, September 2012.
%           CODE paper: M. Landreman et al. Comp. Phys. Comm. 185, 3 (2014)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%this ugly thing handles time-variable input, 
%allowing only one of the input parameters to vary in time, 
%automatically rescaling the other vectors to be of the same size.
%Rewrite it to read
%[refresh_times,ms0,Zs0,Ts0,rhos0,nes0,EHat0]=rescaleParameters(params,refresh_times);
%at some point.
refresh_times = params.refresh_times; 
ms0   = rescaleParameters(params.ms,   refresh_times);
Zs0   = rescaleParameters(params.Zs,   refresh_times);
Ts0   = rescaleParameters(params.Ts,   refresh_times);
rhos0 = rescaleParameters(params.rhos, refresh_times);
nes0  = rescaleParameters(params.nes,  refresh_times);
EHat0 = rescaleParameters(params.EHat, refresh_times);


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
% gridMode             ; 0: Uniform grid with Ny and yMax as given
%                   'auto': Uniform grid with Ny and yMax automatically 
%                           calculated based on tMax to yield well-
%                           converged solution. Does not affect Nt.


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
Ny    = grid.Ny;
Nxi   = grid.Nxi;
yMax  = grid.yMax;


% Boundary condition at yMax.
yMaxBoundaryCondition = 4;
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

%print a couple of relevant physical parameters
fprintf('EHat: %g, Zeff: %g, nbar: %g \n',EHat0(1), sum(rhos0(1,:).*Zs0(1,:)), ms0(1,1)*sum(rhos0(1,:).*Zs0(1,:)./ms0(1,:)) )

startTime = tic;

% Generate differentiation matrices.
% The yWeights vector could be multiplied by any vector of function
% values on the y grid to integrate in y, but we will not need this
% feature for the present application.

yMin=0;
scheme = 12; %12: uniform grid
switch settings.gridMode
    case 0
        [x,~, ddy, d2dy2] = m20121125_04_DifferentiationMatricesForUniformGrid(Ny, yMin, yMax, scheme);
    case 1 %avoid using this setting, experimental non-uniform grid
        scheme = 12; %12: uniform grid
        [s,~, dds, d2ds2] = m20121125_04_DifferentiationMatricesForUniformGrid(Ny, yMin, 1, scheme);
        
        c0 = 1;
        c1 = 0.01;
        c2 = 5;
        x = c1*s.^c2 + c0*s;
        dyds = c2*c1*s.^(c2-1) + c0;
        d2yds2 = c2*(c2-1)*c1*s.^(c2-2);
        
        ddy = diag(1./dyds)*dds;
        d2dy2 = -diag(d2yds2 ./ (dyds.^3)) * dds + diag((1./dyds).^2)*d2ds2;
        d2dy2(end,:) = 0;
    case 'auto' %overrides dx, yMax and Ny to set ''suitable'' values based on the values 
                %chosen for tMax and EHat. Solutions will be well converged (although 
                %it does not set Nxi, which has to be chosen suitably by the user)
        [Ec,xc1,xc2] = runaway_parameters(params);
        yMax = xc2+6;
        dx0 = 0.45;
        dx = dx0/tMax^(1/4);
        Ny = round(yMax/dx);
        if Ny > 1500
            Ny = 1500;
        end
        [x,~, ddy, d2dy2] = m20121125_04_DifferentiationMatricesForUniformGrid(Ny, yMin, yMax, scheme);
end
%print grid parameters
fprintf('Ny: %d,   yMax: %g,   Nxi: %d,   dt: %g,  tMax: %g, Nt: %d\n',Ny,yMax,Nxi,dt,tMax,Nt)

% Make x a row vector:
x = x';


% Order of rows in the matrix and right-hand side:
% --------------------
% for L=0:(Nxi-1)
%   for iy = 1:(Ny-2)
%     Impose kinetic equation
%   Impose boundary condition at yMax
% Enforce regularity: dF/dy=0 at y=0

% Order of columns in the matrix, corresponding to rows in the solution vector:
% --------------------
% for L=0:(Nxi-1)
%   for iy = 1:(Ny-1)
%     Value of F
% Value of F at y=0

matrixSize = Nxi*(Ny-1) + 1;

% Predict roughly how many nonzero elements will be in the sparse
% matrix. This speeds up the code by eliminating the need to reallocate
% memory during matrix construction.
% It is sensitive to which self-collision operator is used -- the 
% conserving ones add dense blocks.
predictedFillFactor = (3*Nxi*nnz(abs(ddy) + abs(d2dy2)))/(matrixSize*matrixSize);
elseif settings.momentumConservation && settings.energyConservation
    predictedFillFactor = predictedFillFactor + 2*Ny*Ny/(matrixSize*matrixSize);
elseif settings.momentumConservation || settings.energyConservation
    predictedFillFactor = predictedFillFactor + Ny*Ny/(matrixSize*matrixSize);
end

fprintf('Matrix size: %g\n',matrixSize)

f = zeros(matrixSize, Nt);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial condition:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fMinus1 = zeros(matrixSize,1);
switch settings.initialDistribution
    case 0
        % A stationary Maxwellian. This sets our normalization of f
        fMinus1(1:(Ny-2)) = exp(-x(2:(Ny-1)).^2);
        fMinus1(end) = 1;
    case 1
        % A shifted (non-normalized) Maxwellian
        x0 = .3; %flow velocity
        fMinus1(1:Ny-2) = (1-x0^2).*exp(-x(2:Ny-1).^2);
        fMinus1(Ny:2*Ny-3) = 2*x0 * x(2:Ny-1).*exp(-x(2:Ny-1).^2);
        fMinus1(end) = 1;
    otherwise
        error('Invalid setting for initial distribution')
end  
fMinus2 = fMinus1;
f(:,1) = fMinus1;

ne0 = nes0(1,1);
Ta0 = Ts0(1,1);
ma0 = ms0(1,1);
me = 1/1822.89;
kappa_e = sqrt(me*Ta0/ma0);

Phi = @(r) erf(r); %Error function = 2/sqrt(pi) int_0^r ds exp(-s^2)
dPhi = @(r) 2/sqrt(pi) * exp(-r.^2); 
G = @(r)(Phi(r)-r.*dPhi(r))./(2*r.^2); %Chandrasekhar function
dG =@(r) 2/sqrt(pi) * (1+1./(r.^2)).*exp(-r.^2) - Phi(r)./(r.^3);

x2=x.*x;

refresh_counter = 0;
timeBeforeTimeAdvance = tic;
for iteration = 2:Nt
    if length(refresh_times) > refresh_counter  %check if it's time to 
                                                %refresh the matrix.
        if refresh_times(refresh_counter + 1) <= tHat(iteration)
            refresh_counter = refresh_counter + 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Begin building matrix.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ms = ms0(refresh_counter,:);
            Ts = Ts0(refresh_counter,:);
            rhos = rhos0(refresh_counter,:);
            Zs = Zs0(refresh_counter,:);
            ne = nes0(refresh_counter);
            EHat = EHat0(refresh_counter);
            
            kappas = sqrt(ms*Ta0./(ma0*Ts));
            lnLambda0 = 1;
            lnLambda = 1;
            collision_time_mod = ne * lnLambda / (ne0 * lnLambda0);

            %I'm sorry that this is a bit ugly -- lots of if statements in
            %this part of the code due to various settings, could possibly
            %be compactified
            if settings.electronCollisions
                coefficientOf_f = 2*Ta0*(2./x.*G(kappa_e*x) + kappa_e*dG(kappa_e*x));
                coefficientOf_dfdy = 2*Ta0*G(kappa_e*x) + G(kappa_e*x)./x2 + kappa_e*dG(kappa_e*x)./x;
                coefficientOf_d2fdy2 = G(kappa_e*x)./x;
            else
                coefficientOf_f = 0;
                coefficientOf_dfdy = 0;
                coefficientOf_d2fdy2 = 0;
            end
               
            
            if settings.approximateCollOp
                coefficientOf_f = zeros(size(x));
                coefficientOf_dfdy = (nbar./x.^2 - nbar./(2*x.^4));
                coefficientOf_d2fdy2 = nbar./(2*x.^3);
            else
                for i=1:length(rhos)
                    coefficientOf_f = coefficientOf_f + rhos(i)*Zs(i)*2*Ta0/Ts(i)*(2*G(kappas(i)*x)./x + kappas(i)*dG(kappas(i)*x));
                    coefficientOf_dfdy = coefficientOf_dfdy + rhos(i)*Zs(i)*(2*Ta0/Ts(i)*G(kappas(i)*x) + G(kappas(i)*x)./x2 + kappas(i)*dG(kappas(i)*x)./x);
                    coefficientOf_d2fdy2 = coefficientOf_d2fdy2 + rhos(i)*Zs(i)*G(kappas(i)*x)./x;
                end
            end
           

            energyScattering = (diag(coefficientOf_f) ...
                + diag(coefficientOf_d2fdy2)*d2dy2 + diag(coefficientOf_dfdy)*ddy);
            if yMaxBoundaryCondition==4
                fakeViscosity = exp((x-yMax)/(0.1));
                energyScattering = energyScattering + (1e-2)*diag(fakeViscosity)*d2dy2;
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
            energyConservingLittleMatrix = zeros(Ny,Ny);
            if settings.energyConservation
                energyConservingLittleMatrix = generateEnergyConservingLittleMatrix(Ta0,Ts,Phi,dPhi,collision_time_mod,rhos,Zs,G,x);
            end
            momentumConservingLittleMatrix = zeros(Ny,Ny);
             if settings.momentumConservation
                momentumConservingLittleMatrix = generateMomentumConservingLittleMatrix(Ta0,Ts,Phi,dPhi,collision_time_mod,rhos,Zs,G,x);
            end


            % Initialize arrays for building the sparse matrix:
            sparseCreatorIndex=1;
            estimated_nnz = 0;
            sparseCreator_i=0;
            sparseCreator_j=0;
            sparseCreator_s=0;
            resetSparseCreator(matrixSize,predictedFillFactor)

            if yMaxBoundaryCondition==3
                rowRange = 1:(Ny-1);
            else
                rowRange = 1:(Ny-2);
            end

            for L=0:(Nxi-1)
                % Add collision operator
                xPartMatrix = energyScattering - L*(L+1) * xPartOfPitchAngleScatteringMatrix;
                if L==0
                    xPartMatrix = xPartMatrix + energyConservingLittleMatrix;
                end
                if L==1
                    xPartMatrix = xPartMatrix + momentumConservingLittleMatrix;
                end

                rowIndices = L*(Ny-1) + rowRange;
                columnIndices = L*(Ny-1) + (1:(Ny-1));
                addSparseBlock(rowIndices, columnIndices, xPartMatrix(1+rowRange, 2:(Ny)))
                if L==0
                    addSparseBlock(rowIndices, matrixSize, xPartMatrix(1+rowRange, 1))
                end


                % Add electric field term- 
                % Sub-diagonal term in L:
                ell = L-1;
                if ell >=0
                    columnIndices = ell*(Ny-1) + (1:(Ny-1));
                    littleMatrix = L/(2*L-1)* EHat*ddy  - diag(EHat*(L-1)*L/(2*L-1)./x);
                    addSparseBlock(rowIndices, columnIndices, littleMatrix(1+rowRange, 2:(Ny)))
                    if ell==0
                        addSparseBlock(rowIndices, matrixSize, littleMatrix(1+rowRange, 1))
                    end
                end

                % Add electric field term- 
                % Super-diagonal term in L:
                ell = L+1;
                if ell<Nxi
                    columnIndices = ell*(Ny-1) + (1:(Ny-1));
                    littleMatrix = (L+1)/(2*L+3)*EHat*ddy  + diag(EHat*(L+1)*(L+2)/(2*L+3)./x);
                    addSparseBlock(rowIndices, columnIndices, littleMatrix(1+rowRange, 2:(Ny)))
                end

            end

            operator = createSparse();
            resetSparseCreator(matrixSize,predictedFillFactor)

            indices = 1:(matrixSize-1);
            addToSparse(indices, indices, ones(size(indices)))

            % For the special point at y=0, apply Neumann condition dF/dy=0:
            addToSparse(matrixSize, matrixSize, ddy(1,1))
            addSparseBlock(matrixSize, 1:(Ny-1), ddy(1, 2:Ny))

            % Impose boundary conditions at y=yMax:
            if yMaxBoundaryCondition==2
                % Add Robin boundary condition dF/dy + (2/y)*F = 0 at yMax:
                for L=0:(Nxi-1)
                    rowIndex = L*(Ny-1) + Ny-1;
                    columnIndices = L*(Ny-1) + (1:(Ny-1));
                    boundaryCondition = ddy(Ny,:);
                    boundaryCondition(Ny) = boundaryCondition(Ny) + 2/yMax - 1; %Subtract 1 since we already put a 1 on the diagonal.
                    addSparseBlock(rowIndex, columnIndices, boundaryCondition(2:Ny))
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

            timeToAssembleMatrix = toc(startTime);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % End of building the matrix.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tic
            [factor_L, factor_U, factor_P, factor_Q] = lu(timeAdvanceMatrix);
            timeToLUFactorize = toc;

    
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

    % Handle boundary condition at y=0:
    rhs(end) = 0;

    % Handle boundary condition at y = yMax:
    if yMaxBoundaryCondition ~= 3
        rhs((1:Nxi)*(Ny-1)) = 0;
    end

    % Now step forward in time.
    % The next line is equivalent to 'soln = matrix \ rhs', but much faster:
    soln = factor_Q * (factor_U \ (factor_L \ (factor_P * rhs)));
    f(:,iteration) = soln;

    fMinus2 = fMinus1;
    fMinus1 = soln;
    
end

fprintf('Done.\n')
fprintf('Time for matrix assembly: %gs,  LU factorization: %gs,  time-advance: %gs.\n',timeToAssembleMatrix, timeToLUFactorize, toc(timeBeforeTimeAdvance))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of time-advance loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The functions below are all utilities for building sparse matrices:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: make the interface clear: which variables are written two
function resetSparseCreator(matrixSize,predictedFillFactor)
    sparseCreatorIndex=1;
    estimated_nnz = floor(matrixSize*matrixSize*predictedFillFactor);
    sparseCreator_i=zeros(estimated_nnz,1);
    sparseCreator_j=zeros(estimated_nnz,1);
    sparseCreator_s=zeros(estimated_nnz,1);
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
    resetSparseCreator()
end


function A = rescaleParameters(X,refresh_times)
    N_rfrs = length(refresh_times);
    [N_rows, N_cols] = size(X);
    if N_rows < N_rfrs
        %fprint('Size of X not compatible with number of refresh times. Setting X constant...')
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
function MC_little_matrix = generateMomentumConservingLittleMatrix(Ta0, Ts, Phi, dPhi, collision_time_mod, rhos, Zs, G, x)
    fM = exp( - x'.^2 * Ta0/Ts(1));
    ws = simpson_quad(x);
    U_DOWN = ws .* x .* ( Phi(x) - x.*dPhi(x) ) * fM;
    Uterm_cols = ws .* ( Phi(x) - x.*dPhi(x) ) / U_DOWN;
    Uterm_rows = collision_time_mod * rhos(1) * Zs(1) * 4*G(x) .* fM';

    MC_little_matrix = Uterm_rows'*Uterm_cols;
end
    
% TODO: package the input variables into one record to reduce the number of arguments
function EC_little_matrix = generateEnergyConservingLittleMatrix(Ta0, Ts, Phi, dPhi, collision_time_mod, rhos, Zs, G, x)
    fM = exp( - x'.^2 * Ta0/Ts(1));
    ws = simpson_quad(x);
    Q_DOWN = ws .* (2*x.^3) .* ( Phi(x) - 2*x.*dPhi(x) ) * fM;
    Qterm_cols = ws .* (2*x) .* ( Phi(x) - 2*x.*dPhi(x) ) / Q_DOWN;
    Qterm_rows = collision_time_mod * rhos(1) * Zs(1) * 2./x .*(Phi(x) - 2*x.*dPhi(x)).* fM';
    Qterm_rows(1) = collision_time_mod * rhos(1) * Zs(1) * (-4/sqrt(pi));
    
    EC_little_matrix = Qterm_rows'*Qterm_cols;
end

end

