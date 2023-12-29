function ns = numOfBasis(c, n)
%
% c: chosen order of Legendre basis
% n: number of dimensions
% 
% for 2D case: ns = (c + 1)*(c + 2)/2; 
%
    ns = factorial(n + c) / (factorial(c) * factorial(n));
end
