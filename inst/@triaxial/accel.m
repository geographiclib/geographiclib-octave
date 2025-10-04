function [v, acc, K] = accel(t, r, v)
%ACCEL  Compute the acceleration of an object constrained to the ellipsoid
%
%   [V, acc, K] = ACCEL(t, R, V)
%
%   Input:
%     t the triaxial ellipsoid object
%     R an n x 3 array of cartesian points on the ellipsoid
%     V an n x 3 array of velocities
%   Output:
%     V the input velocity corrected to be tangent to the ellipsoid
%     acc an n x 3 array of accelerations (normal to the surface)
%     K the Gaussian curvature
%
%   This routine is used for the integration of the geodesic equations.  R
%   should lie on the surface and V should be a unit vector tangent to the
%   surface.  However this function enforces this by an initial call to
%   CART2NORM in order to make the solution of geodesic ODEs better behaved.
%   K is returned in order to allow the ODEs for the reduced length and
%   geodesic scaled to be solved.
%
%   See also CART2NORM, HYBRID, RECKON, DISTANCE

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

% The expression of the acceleration can be found in Panou and Korakitis,
% J. Geod. Sci. (2019).  The expression for K is from
% https://www.johndcook.com/blog/2019/10/07/curvature-of-an-ellipsoid/
  [r, v] = cart2norm(t, r, v);
  axes2 = t.axes.^2;
  n = size(r, 1);
  if n == 1
    up = r ./ axes2;
    u2 = up * up';
    acc = -((v.^2 * (1./axes2)') / u2) * up;
    K = 1 / (prod(axes2) * u2^2);
  else
    up = r ./ axes2;                    % n x 3
    u2 = vecdot(up, up);                % n x 1
    acc = -vecdot(v.^2, 1./axes2) ./ u2 .* up;
    K = 1 ./ (prod(axes2) * u2.^2);
  end
end
