function geod = elliptogeod(t, ellip)
%ELLIPTOGEOD  Convert a general point from ellipsoidal to geodetic
%
%   geod = ELLIPTOGEOD(t, ellip)
%   geod3 = ELLIPTOGEOD(t, ellip3)
%
%   Input:
%     t the trixial ellipsoid object
%     ellip an n x 2 array of ellipsoidal coordinates [bet, omg]
%     ellip3 an n x 3 array of ellipsoidal coordinates [bet, omg, u]
%   Output:
%     geod an n x 2 array of geodetic coordinates [phi, lam]
%     geod3 an n x 3 array of geodetic coordinates [phi, lam, h]
%
%   bet, omg, phi, and lam are measured in degrees.  This routine calls
%   ELLIPTOCART and either CARTTOGEOD or CART2TOGEOD
%
%   See also GEODTOELLIP, ELLIPTOCART, ELLIPTOCART2, CARTTOGEOD, CART2TOGEOD

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = elliptocart(t, ellip);
  if size(ellip, 2) == 3
    geod = carttogeod(t, r);
  else
    geod = cart2togeod(t, r);
  end
end
