function [r, v] = cart2rand(t, n)
%CART2RAND  Return random points on the ellipsoid.
%
%   R = CART2RAND(t, n)
%   [R, V] = CART2RAND(t, n)
%
%   Input:
%     t the triaxial ellipsoid object
%     n the number of points needed
%   Output:
%     R an n x 3 array of cartesian points
%     V an n x 3 array of cartesian velocities
%
%   R are uniformly distributed on the surface of the ellipsoid.  V are unit
%   vectors uniformly distributed in the plane tangent to the surface at R.
%
%   See also CART2NORM

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

% This uses the simple rejection technique given by Marples and Williams,
% Num. Alg. (2023), Algorithm 1 based on the general method of Williamson,
% Phys. Med. Biol. (1987).
  if nargin < 2, n = 1; end
  rn = vecunit(randn(n, 3));
  r = rn .* t.axes;
  up = rn ./ t.axes;
  g = t.c * vecabs(up);
  rej = rand(n, 1) > g;
  if any(rej)
    r(rej,:) = cart2rand(t, sum(rej));
  end
  if nargout > 1
    v = randn(n, 3);
    % up(rej,:) = r(rej,:) ./ t.axes.^2
    up = r ./ t.axes.^2;
    u2 = sum(up.^2, 2);
    v = vecunit(v - up .* vecdot(up, v) ./ u2);
  end
end
