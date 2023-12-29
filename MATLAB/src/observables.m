function H = observables(ns,Nx,ind,MLPC)
%
% Observables
%
% ns: total number of basis
% Nx: total number of variables
% ind: Legendre basis index
% MLPC: Coefficients of multi-variate Legendre polynomial
%
    H = zeros(Nx,ns);
    par = zeros(1,Nx);

    for i = 1:Nx
        % Observable Polynomials
        Obs = zeros(1,Nx+1);
        Obs(1,1) = 1;
        Obs(1,i+1) = 1;
    
        % Matix integration
        for k = 1:ns
            for j = 1:ns
                if (MLPC(k,j) ~= 0)
                    flag = 1;
                    for in = 1:Nx
                        par(in) = round(ind(j,in) + Obs(1,in+1) + 1);
                        if (mod(par(in),2) == 0)
                            flag = 0;
                            break;
                        end
                    end
                    if (flag == 1)
                        H(i,k) = H(i,k) + (2^Nx)*MLPC(k,j)*Obs(1,1) / (prod(par));
                    end
                end
            end
        end
    end

end