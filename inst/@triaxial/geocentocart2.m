function r = geocentocart2(t, geocen)
%GEOCENTOCART2  Convert a  point from geocentric to cartesian
%
%   r = GEOCENTOCART2(t, geocen)
%
%   Input:
%     t the triaxial ellipsoid object
%     geocen an n x 2 array of the geocentric coordinates [phic, lamc]
%   Output:
%     r an n x 3 array of cartesian points lying on the ellipsoid
%
%   phic and lamc are measured in degrees.
%
%   See also CART2TOGEOCEN

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = [cosd(geocen(:,1)) .* [cosd(geocen(:,2)), sind(geocen(:,2))], ...
       sind(geocen(:, 1))];
  r = r ./ vecabs(r./t.axes);
end
