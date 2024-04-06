function param = geodtoparam(t, geod)
%GEODTOPARAM  Convert a general point from parametric to geodetic
%
%   param = GEODTOPARAM(t, geod)
%   param3 = GEODTOPARAM(t, geod3)
%
%   Input:
%     t the trixial ellipsoid object
%     geod an n x 2 array of geodetic coordinates [phi, lam]
%     geod3 an n x 3 array of geodetic coordinates [phi, lam, h]
%   Output:
%     param an n x 2 array of parametric coordinates [phip, lamp]
%     param3 an n x 3 array of parametric coordinates [phip, lamp, h]
%
%   phi, lam, phio, and lamp are measured in degrees.  This routine calls
%   GEODTOCART and CART2TOPARAM.
%
%   See also PARAMTOGEOD, GEODTOCART, CART2TOPARAM, CARTTOPARAM

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = geodtocart(t, geod(:, 1:2));
  param = cart2toparam(t, r);
  if size(geod, 2) == 3
    param = [param, geod(:, 3)];
  end
end
