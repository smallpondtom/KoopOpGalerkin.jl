function ind = basisIndex(c, ns, Nx)
% 
% c: chosen order of polynomial
% ns: total number of basis functions
% Nx: total number of variables 
%
    ind = zeros(ns, Nx);
    s = 1;
    for ord = 0:c % Start from 0 to include the 0th order
        [ind, s] = iterateDims(ord, zeros(1, Nx), 1, c, Nx, ind, s);
    end
end

function [ind, s] = iterateDims(ord, currentComb, currentDim, c, Nx, ind, s)
%
% This implementation uses a helper function `iterateDims` to handle the 
% recursion. The function `iterateDims` iterates over the indices for each 
% dimension and calls itself for the next dimension. It stops when it 
% reaches the last dimension and updates the ind matrix with the current 
% combination of indices.
%
    if currentDim == Nx
        currentComb(Nx) = ord - sum(currentComb(1:Nx-1));
        if sum(currentComb) == ord
            ind(s, :) = currentComb;
            s = s + 1;
        end
        return;
    end

    for i = 0:ord
        newComb = currentComb;
        newComb(currentDim) = i;
        if sum(newComb) <= ord
            [ind, s] = iterateDims(ord, newComb, currentDim + 1, c, Nx, ind, s);
        end
    end
end

%% Version 2 - Recursive
% function ind = basisIndex(c, ns, Nx)
%     ind = zeros(ns, Nx);
%     [ind, ~] = iterateDims(1, zeros(1, Nx), 1, c, Nx, ind);
% end
% 
% function [ind, s] = iterateDims(s, currentComb, currentDim, c, Nx, ind)
%     if currentDim > Nx
%         % Check if the sum of the current combination is within the order.
%         if sum(currentComb) <= c
%             ind(s, :) = currentComb;
%             s = s + 1;
%         end
%         return;
%     else
%         for i = 0:c
%             newComb = currentComb;
%             newComb(currentDim) = i;
%             [ind, s] = iterateDims(s, newComb, currentDim + 1, c, Nx, ind);
%         end
%     end
% end


%% Version 3 - Recursive (different order)
% function ind = basisIndex(c, ns, Nx)
%     ind = zeros(ns, Nx);
%     [ind, ~] = iterateOrder(1, 1, c, Nx, ind);
% end
% 
% function [ind, s] = iterateOrder(s, currentOrder, c, Nx, ind)
%     if currentOrder > c
%         return;
%     end
%     for i = 0:currentOrder
%         currentComb = zeros(1, Nx);
%         currentComb(1) = i;
%         if Nx > 1
%             currentComb(2) = currentOrder - i;
%         end
%         if sum(currentComb) == currentOrder
%             ind(s, :) = currentComb;
%             s = s + 1;
%         end
%         if Nx > 2
%             [ind, s] = iterateDims(s, currentComb, 3, currentOrder - sum(currentComb), c, Nx, ind);
%         end
%     end
%     [ind, s] = iterateOrder(s, currentOrder + 1, c, Nx, ind);
% end
% 
% function [ind, s] = iterateDims(s, currentComb, currentDim, remainingOrder, c, Nx, ind)
%     if currentDim > Nx
%         return;
%     end
%     for i = 0:min(remainingOrder, c)
%         newComb = currentComb;
%         newComb(currentDim) = i;
%         if sum(newComb) <= c
%             if currentDim == Nx
%                 ind(s, :) = newComb;
%                 s = s + 1;
%             else
%                 [ind, s] = iterateDims(s, newComb, currentDim + 1, remainingOrder - i, c, Nx, ind);
%             end
%         end
%     end
% end
