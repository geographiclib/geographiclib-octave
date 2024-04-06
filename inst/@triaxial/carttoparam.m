function param3 = carttoparam(t, r)
%CARTTOPARAM  Convert a surface point from cartesion to parametric
%
%   param3 = CARTTOPARAM(t, r)
%
%   Input:
%     t the trixial ellipsoid object
%     r an n x 3 array of cartesian points
%   Output:
%     param3 an n x 3 array of parametric coordinates [phip, lamp, h]
%
%   phip and lamp are measured in degrees.  This routine calls CARTTOCART2
%   follwoed by CART2TOPARAM.
%
%   See also CARTTOCART2, CART2TOPARAM, PARAMTOCART2, PARAMTOCART

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  [r, h] = carttocart2(t, r);
  param = cart2toparam(t, r);
  param3 = [param, h];
end
