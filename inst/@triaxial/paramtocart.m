function r = paramtocart(t, param)
%PARAMTOCART  Convert a general point from parametric to cartesian
%
%   r = PARAMTOCART(t, param)
%   r = PARAMTOCART(t, param3)
%
%   Input:
%     t the trixial ellipsoid object
%     param an n x 2 array of the parametric coordinates [phip, lamp]
%     param3 an n x 3 array of the parametric coordinates [phip, lamp, h]
%   Output:
%     r an n x 3 array of cartesian points
%
%   phip and lamp are measured in degrees.  With param (an n x 2 array), h is
%   assumed to be 0 so that r lies on the ellipsoid.
%
%   See also CARTTOPARAM, CART2TOCART

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  cp = cosd(param(:,1));
  r = t.axes .* ...
      [cp .* cosd(param(:,2)), cp .* sind(param(:,2)), sind(param(:, 1))];
  if size(param, 2) == 3
    r = cart2tocart(t, r, param(:,3));
  end
  % North and east vectors proportional to dx/dphi and dx/dlam
  % N:[a,b,c]*[-sin(phip)*cos(lamp),-sin(phip)*sin(lamp),cos(phip)];
  % E:[a,b,c]*[-cos(phip)*sin(lamp), cos(phip)*cos(lamp), 0];
  % N.B. these are NOT orthogonal, N.E =
  % NE:(a^2-b^2) * cos(phip)*sin(phip)*cos(lamp)*sin(lamp);
end
