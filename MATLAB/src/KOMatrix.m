function K = KOMatrix(ns,Nx,Nt,fx,ind,NLPC,DNLPC,MLPC)
%
% Compute Koopman Operator Matrix (KOMatrix)
% GENERALIZED CODE: For Nx variables
%
% ns: total number of basis
% Nx: total number of variables
% Nt: total number of terms in f(x)
% fx: function array including coefficients and exponent info
% ind: Legendre basis index
% NLPC: coefficients of normalized Legendre polynomials
% DNLPC: coefficients of derivative of normalized Legendre Polynomials
% MLPC: coefficients of multi-variate Legendre polynomials
%
    K = zeros(ns,ns);
    MDLP = zeros(Nx,ns);
    par = zeros(1,Nx);
    
    % i represents the selected Legendre polynomials we are evaluating 
    % the derivatives of
    for i = 1:ns
    
        % Partial Derivatives of Legendre polynomials in multiple dimensions
        for j = 1:ns
            % Partial derivatives are chain rule so just multiply 
            % e.g., dL1dt * L2 * L3 * ... * L_{Nx}
            % or,   L1 * dL2dt * L3 * ... * L_{Nx}
            for k = 1:Nx
                MDLP(k,j) = DNLPC(ind(i,k)+1,ind(j,k)+1);
                for l = 1:Nx
                    if k ~= l
                        MDLP(k,j) = MDLP(k,j) * NLPC(ind(i,l)+1,ind(j,l)+1);
                    end
                end
            end
        end
    
        % Total Derivative = fx * grad(Legendre Polynomials)
        DB = zeros(ns*2,Nx+1);
        s = 0;
        for dim = 1:Nx
            for ifx = 1:Nt  % loop thru all the terms in f(x)
                if (fx(dim,ifx,1) == 0)
                    break;
                else
                    for j = 1:ns
                        if (MDLP(dim,j) ~= 0)
                            s = s + 1;
                            DB(s,1) = fx(dim,ifx,1)*MDLP(dim,j);
                            % For each variable update the database
                            for k = 2:Nx+1
                                DB(s,k) = fx(dim,ifx,k) + ind(j,k-1);
                            end
                        end
                    end
                end
            end
        end
    
        % Matix integration
        for k = 1:ns
            for j = 1:ns
                if (MLPC(k,j) ~= 0)
                    for ifx = 1:s
                        flag = 1;
                        for in = 1:Nx
                            % First part of integration where we get the
                            % (exponent+1) with which we divide the integrand
                            par(in) = round(ind(j,in) + DB(ifx,in+1) + 1);
                            % Check if exponent is "even" or "odd"
                            if (mod(par(in),2) == 0)
                                flag = 0;
                                break;
                            end
                        end
                        % If odd, perform the second part of the integration
                        if (flag == 1)
                            K(i,k) = K(i,k) + (2^Nx)*MLPC(k,j)*DB(ifx,1) / (prod(par));
                        end
                    end
                end
            end
        end
    end

end