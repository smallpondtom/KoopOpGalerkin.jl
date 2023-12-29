function LPC = LP(c)
%
% Coefficients (C) of the Legendre Polynomials (LP)
% c: order of polynomial 0, 1, 2, ..., c
%
    % so dimension becomes (c+1)
    LPC = zeros(c+1,c+1);
    LPC(1,1) = 1;
    LPC(2,2) = 1;
    for i = 3:c+1
        for j = 1:i-1
            LPC(i,j+1) = LPC(i,j+1) + (2*(i-2) + 1)/(i-1)*LPC(i-1,j);
            LPC(i,j) = LPC(i,j) - (i-2)/(i-1)*LPC(i-2,j);
        end
    end
end