function r = paramtocart2(t, param)
%PARAMTOCART2  Convert a point from parametric to cartesian
%
%   r = PARAMTOCART2(t, param)
%
%   Input:
%     t the triaxial ellipsoid object
%     param an n x 2 array of the parametric coordinates [phip, lamp]
%   Output:
%     r an n x 3 array of cartesian points lying on the ellipsoid
%
%   phip and lamp are measured in degrees.
%
%   See also CART2TOPARAM

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = t.axes .* ...
      [cosd(param(:,1)) .* [cosd(param(:,2)), sind(param(:,2))], ...
       sind(param(:, 1))];
  % r = r ./ vecabs(r./t.axes)   ...  a no-op
  % North and east vectors proportional to dx/dphi and dx/dlam
  % N:[a,b,c]*[-sin(phip)*cos(lamp),-sin(phip)*sin(lamp),cos(phip)];
  % E:[a,b,c]*[-cos(phip)*sin(lamp), cos(phip)*cos(lamp), 0];
  % N.B. these are NOT orthogonal, N.E =
  % NE:(a^2-b^2) * cos(phip)*sin(phip)*cos(lamp)*sin(lamp);
end
