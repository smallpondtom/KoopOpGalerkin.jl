function NLPC = NLP(c,LPC)
%
% Coefficients (C) of the Normalized Legendre Polynomials (NLP)
%
% c: chosen order of the polynomial
% LPC: Coefficients of Legendre polynomial
%
    NLPC = zeros(c+1,c+1);
    for i = 1:c+1
        for j = 1:c+1
            NLPC(i,j) = sqrt((2*(i-1)+1)/2)*LPC(i,j);
        end
    end
end