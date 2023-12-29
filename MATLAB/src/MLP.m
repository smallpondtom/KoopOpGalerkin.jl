function MLPC = MLP(ns,Nx,ind,NLPC)
%
% Coefficients (C) of the Multi-variate Legendre polynomials (MLP)
% ns: total number of basis
% Nx: total number of variables
% ind: Legendre Polynomial index
% NLPC: coefficients of normalized Legendre polynomial
%
    MLPC = ones(ns,ns);
    for i = 1:ns
        for j = 1:ns
            for k = 1:Nx
                MLPC(i,j) = MLPC(i,j) * NLPC(ind(i,k)+1,ind(j,k)+1);
            end
        end
    end
end