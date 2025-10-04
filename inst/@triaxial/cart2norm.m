function [r, v] = cart2norm(t, r, v)
%CART2NORM  Force a point to lie on the surface point of the ellipsoid
%
%   R = CART2NORM(t, R)
%   [R, V] = CART2NORM(t, R, V)
%
%   Input:
%     t the triaxial ellipsoid object
%     R an n x 3 array of cartesian points on the ellipsoid
%     V an n x 3 array of cartesian velocities
%   Output:
%     R an n x 3 array of corrected cartesian points
%     V an n x 3 array of corrected cartesian velocities
%
%   The expectation is that R and V start close to the ellipsoid and this
%   routine merely tweaks these to force R to be on the surface and V to be a
%   unit vector tangent to the surface at R.
%
%   See also CARTTOCART2, CART2TOELLIP, CART2TOGEOD, CART2TOPARAM,
%     CARTTOGEOD, CARTTOPARAM

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  rn = r ./ t.axes;
  ra = vecabs(rn);
  rnx = rn ./ ra;
  % Don't normalize if ra is within 1 +/- eps
  l = abs(ra - 1) > eps;
  rn(l,:) = rnx(l,:);
  r(l,:) = rn(l,:) .* t.axes;
  if nargin < 3, return; end
  up = rn ./ t.axes;
  u2 = sum(up.^2, 2);
  v = vecunit(v - up .* vecdot(up, v) ./ u2);
end
