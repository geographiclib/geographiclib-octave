function param = elliptoparam(t, ellip)
%ELLIPTOPARAM  Convert a general point from ellipsoidal to geodetic
%
%   param = ELLIPTOPARAM(t, ellip)
%   param3 = ELLIPTOPARAM(t, ellip3)
%
%   Input:
%     t the trixial ellipsoid object
%     ellip an n x 2 array of ellipsoidal coordinates [bet, omg]
%     ellip3 an n x 3 array of ellipsoidal coordinates [bet, omg, u]
%   Output:
%     param an n x 2 array of parametic coordinates [phip, lamp]
%     param3 an n x 3 array of parametic coordinates [phip, lamp, h]
%
%   bet, omg, phip, and lamp are measured in degrees.  This routine calls
%   ELLIPTOCART and either CARTTOPARAM or CART2TOPARAM
%
%   See also PARAMTOELLIP, ELLIPTOCART, ELLIPTOCART2, CARTTOPARAM, CART2TOPARAM

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = elliptocart(t, ellip);
  if size(ellip, 2) == 3
    param = carttoparam(t, r);
  else
    param = cart2toparam(t, r);
  end
end
