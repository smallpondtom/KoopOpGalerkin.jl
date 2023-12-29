function Sol = KoopOpSol(ord,fx,x0,tspan)
% 
% Solution of nonlinear system using Koopman Operator with Orthogonal
% Polynomial basis.
%
% ord: chosen order of polynomial
% x0: state initial conditions
% tspan: time span of simulation
% fx: function array including coefficients and exponent info
%
    [Nx,Nt,~] = size(fx);
    assert(ord >= 3, 'Order of basis should be greater than 3 for this application.');
    
    %% Number of basis functions
    ns = numOfBasis(ord,Nx);
    
    %% Legendre basis indexes
    ind = basisIndex(ord,ns,Nx);
    
    %% Definition of the Legendre polynomials
    LPC = LP(ord);
    
    %% Legendre polynomials in multiple dimensions
    NLPC = NLP(ord,LPC);
    MLPC = MLP(ns,Nx,ind,NLPC);
    
    %% Derivative of Normalized Legendre polynomials
    DNLPC = DNLP(ord,NLPC);
    
    %% Operator Matrix
    K = KOMatrix(ns,Nx,Nt,fx,ind,NLPC,DNLPC,MLPC);
    
    %% Observables
    H = observables(ns,Nx,ind,MLPC);
    
    %% Eigenvalue decomposition
    [V,D,Wt] = eig(K,'nobalance');  % get left eigenvector
    iV = V \ eye(ns);

    
    %% Initial condition mapped onto Legendre
    L0 = funcSpaceIC(ns,x0,ind,MLPC);
    
    %% Computaion of the solution as a function of time
    nt = length(tspan);
    Sol = zeros(Nx,nt);
    for k = 1:nt
        Sol(:,k) = real(H * V * expm(D * tspan(k)) * iV * L0);
    end
end