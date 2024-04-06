function ellip = geodtoellip(t, geod)
%GEODTOELLIP  Convert a general point from geodetic to ellipsoidal
%
%   ellip = GEODTOELLIP(t, geod)
%   ellip3 = GEODTOELLIP(t, geod3)
%
%   Input:
%     t the trixial ellipsoid object
%     geod an n x 2 array of geodetic coordinates [phi, lam]
%     geod3 an n x 3 array of geodetic coordinates [phi, lam, h]
%   Output:
%     ellip an n x 2 array of ellipsoidal coordinates [bet, omg]
%     ellip3 an n x 3 array of ellipsoidal coordinates [bet, omg, u]
%
%   phi, lam, bet, omg are measured in degrees.  This routine calls
%   GEODTOCART and either CARTTOELLIP or CART2TOELLIP.
%
%   See also ELLIPTOGEOD, GEODTOCART, CARTTOELLIP, CART2TOELLIP

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = geodtocart(t, geod);
  if size(geod, 2) == 3
    ellip = carttoellip(t, r);
  else
    ellip = cart2toellip(t, r);
  end
end
