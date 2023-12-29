function g0 = funcSpaceIC(ns,x0,ind,MLPC)
%
% Generation of the functional space initial conditions
%
% ns: total number of basis
% x0: initial condition of states
% ind: Legendre polynomial index
% MLPC: coefficients of multi-variate Legendre polynomial
%   
    % reshape initial condition into row vector
    [ntmp,mtmp] = size(x0);
    if ntmp > mtmp
        x0 = x0';
    end

    g0 = zeros(ns,1);
    arrExp = @(arr,exp) prod(arr .^ exp);
    for i = 1:ns
        for j = 1:ns
            g0(i) = g0(i) + MLPC(i,j) * arrExp(x0,ind(j,:));
        end
    end
end