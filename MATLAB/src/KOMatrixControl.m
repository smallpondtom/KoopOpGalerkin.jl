function N = KOMatrixControl(ns,Mx,Mt,hx,ind,NLPC,DNLPC,MLPC)
%
% Compute Koopman operator matrix w.r.t the control vector fields
% GENERALIZED CODE: For Nx variables
%
% ns: total number of basis
% Mx: total number of variables of control
% Mt: total number of terms in g(x)
% hx: function array of control including coefficients and exponent info
% ind: Legendre basis index
% NLPC: coefficients of normalized Legendre polynomials
% DNLPC: coefficients of derivative of normalized Legendre Polynomials
% MLPC: coefficients of multi-variate Legendre polynomials
%
    N = zeros(ns,ns);
    MDLP = zeros(Mx,ns);
    par = zeros(1,Mx);
    
    % i represents the selected Legendre polynomials we are evaluating 
    % the derivatives of
    for i = 1:ns
    
        % Partial Derivatives of Legendre polynomials in multiple dimensions
        for j = 1:ns
            % Partial derivatives are chain rule so just multiply 
            % e.g., dL1dt * L2 * L3 * ... * L_{Nx}
            % or,   L1 * dL2dt * L3 * ... * L_{Nx}
            for k = 1:Mx
                MDLP(k,j) = DNLPC(ind(i,k)+1,ind(j,k)+1);
                for l = 1:Mx
                    if k ~= l
                        MDLP(k,j) = MDLP(k,j) * NLPC(ind(i,l)+1,ind(j,l)+1);
                    end
                end
            end
        end
    
        % Total Derivative = hx * grad(Legendre Polynomials)
        DB = zeros(ns*2,Mx+1);
        s = 0;
        for dim = 1:Mx
            for ifx = 1:Mt  % loop thru all the terms in h(x)
                if (hx(dim,ifx,1) == 0)
                    break;
                else
                    for j = 1:ns
                        if (MDLP(dim,j) ~= 0)
                            s = s + 1;
                            DB(s,1) = hx(dim,ifx,1)*MDLP(dim,j);
                            % For each variable update the database
                            for k = 2:Mx+1
                                DB(s,k) = hx(dim,ifx,k) + ind(j,k-1);
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
                        for in = 1:Mx
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
                            N(i,k) = N(i,k) + (2^Mx)*MLPC(k,j)*DB(ifx,1) / (prod(par));
                        end
                    end
                end
            end
        end
    end

end