function geod = paramtogeod(t, param)
%PARAMTOGEOD  Convert a general point from parametric to geodetic
%
%   geod = PARAMTOGEOD(t, param)
%   geod3 = PARAMTOGEOD(t, param3)
%
%   Input:
%     t the trixial ellipsoid object
%     param an n x 2 array of parametric coordinates [phip, lamp]
%     param3 an n x 3 array of parametric coordinates [phip, lamp, h]
%   Output:
%     geod an n x 2 array of geodetic coordinates [phi, lam]
%     geod3 an n x 3 array of geodetic coordinates [phi, lam, h]
%
%   phip, lamp, phi, and lam are measured in degrees.  This routine calls
%   PARAMTOCART and CART2TOGEOD
%
%   See also GEODTOPARAM, PARAMTOCART, PARAMTOCART2, CARTTOGEOD, CART2TOGEOD

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = paramtocart(t, param(:, 1:2));
  geod = cart2togeod(t, r);
  if size(param, 2) == 3
    geod = [geod, param(:, 3)];
  end
end
