function Sol = KoopOpBilinear(ord,fx,hx,u,x0,tspan)
% 
% Solution of Bilinear Polynomial system using Koopman Operator with 
% Orthogonal Polynomial basis.
%
% system:
% xdot(t) = A(x) + \sum_{i}^m B_i(x)u_i(t)
%
% ord: chosen order of polynomial
% fx: function array including coefficients and exponent info (3D array)
% hx: state function array associated with control input including 
%     coeffcients and exponent info 
%     (cell array containing 3D arrays or 3D array)
% u: control inputs cell array of function
% x0: state initial conditions
% tspan: time span of simulation
%
    [Nx,Nt,~] = size(fx);
    if iscell(hx)
        M = length(hx); % number of control inputs
        [Mx,Mt,~] = size(hx{1});
    else
        M = 1;
        [Mx,Mt,~] = size(hx);
    end
    assert(ord >= 3, 'Order of basis should be greater than 3 for this application.');
    
    %% Number of basis functions
    ns = numOfBasis(ord,Nx)
    
    %% Legendre basis indexes
    ind = basisIndex(ord,ns,Nx);
    
    %% Definition of the Legendre polynomials
    LPC = LP(ord);
    
    %% Legendre polynomials in multiple dimensions
    NLPC = NLP(ord,LPC);
    MLPC = MLP(ns,Nx,ind,NLPC);

    size(NLPC)
    size(MLPC)
    
    %% Derivative of Normalized Legendre polynomials
    DNLPC = DNLP(ord,NLPC);
    size(DNLPC)
    
    %% Operator Matrix
    K = KOMatrix(ns,Nx,Nt,fx,ind,NLPC,DNLPC,MLPC);

    %% Control Matrix
    N = zeros(ns,ns,M);
    for i = 1:M
        N(:,:,i) = KOMatrixControl(ns,Mx,Mt,hx{i},ind,NLPC,DNLPC,MLPC);
    end

    %% Observables
    H = observables(ns,Nx,ind,MLPC);
    
    %% Eigenvalue decomposition
    [V,D] = eig(K,'nobalance');
    iV = V \ eye(ns);
    size(D)
    V
    
    %% Initial condition mapped onto Legendre
    L0 = funcSpaceIC(ns,x0,ind,MLPC);
    
    %% Computaion of the solution as a function of time
    
    function ztildedot = KoopmanSystem(t,ztilde)
        % make sure to keep only the real parts
        ztildedot = real(K*ztilde);
        for j = 1:M
            ztildedot = ztildedot + real(N(:,:,j)*iV*ztilde*u{j}(t));
        end
    end

    [~,Z] = ode45(@(t,z) KoopmanSystem(t,z), tspan, L0);
    
    % Fix the dimensions so that row = states & col = time
    [rowdim,coldim] = size(Z);
    if rowdim > coldim
        Z = Z';
    end

    Sol = H*Z;  % Map them to the observables
end