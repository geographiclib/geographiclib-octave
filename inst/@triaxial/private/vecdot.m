function z = vecdot(x, y)
% Input:
%  x,y = n x 3 vectors
% Output:
%  z = n x 1 of dot products
%
% Like dot(x,y,2) but unlike the builtin routine, this can deal with N . 1
% or 1 . N cases.
  z = sum(x .* y, 2);
end
