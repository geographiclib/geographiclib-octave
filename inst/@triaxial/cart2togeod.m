function geod = cart2togeod(t, r)
%CART2TOGEOD  Convert a surface point from cartesion to geodetic
%
%   geod = CART2TOGEOD(t, r)
%
%   Input:
%     t the triaxial ellipsoid object
%     r an n x 3 array of cartesian points on the ellipsoid
%   Output:
%     geod an n x 2 array of geodetic coordinates [phi, lam]
%
%   phi and lam are measured in degrees.  This routine assumes that r lie on
%   the surface of the ellipsoid and that v is a unit vector tangent to the
%   ellipsoid at r.  To ensure that this is the case, call CARTNORM.  To
%   convert arbitrary points use CARTTOGEOD.
%
%   See also CARTNORM, CARTTOGEOD, GEODTOCART2, GEODTOCART

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = r ./ t.axes.^2;
  geod = [atan2d(r(:, 3), hypot(r(:, 2), r(:, 1))), ...
          atan2d(r(:, 2), r(:, 1))];
end
