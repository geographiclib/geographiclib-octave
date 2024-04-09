function geod3 = carttogeod(t, r)
%CARTTOGEOD  Convert a surface point from cartesion to geodetic
%
%   geod3 = CARTTOGEOD(t, r)
%
%   Input:
%     t the triaxial ellipsoid object
%     r an n x 3 array of cartesian points
%   Output:
%     geod3 an n x 3 array of geodetic coordinates [phi, lam, h]
%
%   phi and lam are measured in degrees.  This routine calls CARTTOCART2
%   follwoed by CART2TOGEOD.
%
%   See also CARTTOCART2, CART2TOGEOD, GEODTOCART2, GEODTOCART

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  [r, h] = carttocart2(t, r);
  geod = cart2togeod(t, r);
  geod3 = [geod, h];
end
