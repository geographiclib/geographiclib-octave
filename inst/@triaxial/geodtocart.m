function r = geodtocart(t, geod)
%GEODTOCART  Convert a general point from geodetic to cartesian
%
%   r = GEODTOCART(t, geod)
%   r = GEODTOCART(t, geod3)
%
%   Input:
%     t the triaxial ellipsoid object
%     geod an n x 2 array of the geodetic coordinates [phi, lam]
%     geod3 an n x 3 array of the geodetic coordinates [phi, lam, h]
%   Output:
%     r an n x 3 array of cartesian points
%
%   phi and lam are measured in degrees.  With geod (an n x 2 array), the h
%   is assumed to be 0 so that r lies on the ellipsoid.
%
%   See also CARTTOGEOD, CART2TOGEOD, GEODTOCART2

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = t.axes.^2 .* ...
      [cosd(geod(:,1)) .* [cosd(geod(:,2)), sind(geod(:,2))], ...
       sind(geod(:, 1))];
  r = r ./ vecabs(r./t.axes);
  if size(geod, 2) == 3
    r = cart2tocart(t, r, geod(:,3));
  end
  % North and east vectors proportional to dx/dphi and dx/dlam
  % N:c^2*[-sin(phi)*cos(lam), -si(phi)*sin(lam), cos(phi)]*
  % [a^2, b^2, a^2 * cos(lam)^2 + b^2 * sin(lam)^2];
  % E:[-a^2*cos(phi)*sin(lam), b^2*cos(phi)*cos(lam),c^2*sin(phi)*cos(phi)]
  % *[b^2*cos(phi)^2+c^2*sin(phi)^2,
  % a^2*cos(phi)^2+c^2*sin(phi)^2,
  % (a^2-b^2)*cos(phi)*sin(lam)*cos(lam)];
  % N.B. these are NOT orthogonal, N.E =
  % NE:((a^2-b^2)*c^2*cos(lam)*sin(lam)*cos(phi)*sin(phi))*
  % ((a^2+b^2)*c^2*sin(phi)^2+
  %  (a^2*b^2+a^2*c^2*cos(lam)^2+b^2*c^2*sin(lam)^2)*cos(phi)^2);
end
