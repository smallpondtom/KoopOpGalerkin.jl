function DNLPC = DNLP(c,NLPC)
%
% Coefficients (C) of the Derivative of the Normalized Legendre Polynomials
% (DNLP)
%
% c: chosen order of the polynomial
% NLPC: coefficients of normalized Legendre polynomial
%
    DNLPC = zeros(c+1,c+1);
    for i = 2:c+1
        for j = 1:c
            DNLPC(i,j) = j*NLPC(i,j+1);
        end
    end
end