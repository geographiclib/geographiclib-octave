function geocen = cart2togeocen(~, r)
%CART2TOGEOCEN  Convert a surface point from cartesion to geocentric
%
%   geocen = CART2TOGEOCEN(t, r)
%
%   Input:
%     t the triaxial ellipsoid object
%     r an n x 3 array of cartesian points on the ellipsoid
%   Output:
%     geocen an n x 2 array of geocenraphic coordinates [phic, lamc]
%
%   phic and lamc are measured in degrees.  This routine assumes that r lie on
%   the surface of the ellipsoid and that v is a unit vector tangent to the
%   ellipsoid at r.  To ensure that this is the case, call CARTNORM.  To
%   convert arbitrary points use CARTTOGEOCEN.
%
%   See also CARTNORM, GEOCENTOCART2

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  geocen = [atan2d(r(:, 3), hypot(r(:, 2), r(:, 1))), ...
            atan2d(r(:, 2), r(:, 1))];
end
