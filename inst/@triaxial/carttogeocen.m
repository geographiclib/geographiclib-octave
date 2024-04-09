function geocen3 = carttogeocen(t, r)
%CARTTOGEOCEN  Convert a surface point from cartesion to geocentric
%
%   geocen3 = CARTTOGEOCEN(t, r)
%
%   Input:
%     t the trixial ellipsoid object
%     r an n x 3 array of cartesian points
%   Output:
%     geocen3 an n x 3 array of geocentric coordinates [phic, lamc, h]
%
%   phic and lamc are measured in degrees.  This routine calls CARTTOCART2
%   follwoed by CART2TOGEOCEN.
%
%   See also CARTTOCART2, CART2TOGEOCEN, GEOCENTOCART2, GEOCENTOCART

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  [r, h] = carttocart2(t, r);
  geocen = cart2togeocen(t, r);
  geocen3 = [geocen, h];
end
