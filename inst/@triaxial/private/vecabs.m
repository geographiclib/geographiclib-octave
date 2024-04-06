function y = vecabs(x)
% Input:
%  x = n x 3 vectors
% Output:
%  y = n x 1 of absolute values
%
% Equivalent to vecnorm(x,2,2)
  y = sqrt(sum(x.^2, 2));
end
