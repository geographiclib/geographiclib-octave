function ellip = paramtoellip(t, param)
%PARAMTOELLIP  Convert a general point from parametric to ellipsoidal
%
%   ellip = PARAMTOELLIP(t, param)
%   ellip3 = PARAMTOELLIP(t, param3)
%
%   Input:
%     t the trixial ellipsoid object
%     param an n x 2 array of parametric coordinates [phip, lamp]
%     param3 an n x 3 array of parametric coordinates [phip, lamp, h]
%   Output:
%     ellip an n x 2 array of ellipsoidal coordinates [bet, omg]
%     ellip3 an n x 3 array of ellipsoidal coordinates [bet, omg, u]
%
%   phip, lamp, bet, omg are measured in degrees.  This routine calls
%   PARAMTOCART and either CARTTOELLIP or CART2TOELLIP.
%
%   See also ELLIPTOPARAM, PARAMTOCART, CARTTOELLIP, CART2TOELLIP

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = paramtocart(t, param);
  if size(param, 2) == 3
    ellip = carttoellip(t, r);
  else
    ellip = cart2toellip(t, r);
  end
end
