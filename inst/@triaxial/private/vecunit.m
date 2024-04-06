function y = vecunit(x)
% Input:
%  x = n x 3 vectors
% Output:
%  y = n x 3 vectors converted to unit vectors
  z = vecabs(x);
  z(z == 0) = 1;
  y = x ./ z;
end
