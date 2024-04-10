function param = cart2toparam(t, r)
%CART2TOPARAM  Convert a surface point from cartesion to parametric
%
%   param = CART2TOPARAM(t, r)
%
%   Input:
%     t the triaxial ellipsoid object
%     r an n x 3 array of cartesian points on the ellipsoid
%   Output:
%     param an n x 2 array of parametric coordinates [phip, lamp]
%
%   phip and lamp are measured in degrees.  This routine assumes that r lie on
%   the surface of the ellipsoid and that v is a unit vector tangent to the
%   ellipsoid at r.  To ensure that this is the case, call CARTNORM.  To
%   convert arbitrary points use CARTTOPARAM.
%
%   See also CARTNORM, PARAMTOCART2

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = r ./ t.axes;
  param = [atan2d(r(:, 3), hypot(r(:, 2), r(:, 1))), ...
           atan2d(r(:, 2), r(:, 1))];
end
