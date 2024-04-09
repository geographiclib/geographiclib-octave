function r = geocentocart(t, geocen)
%GEOCENTOCART  Convert a general point from geocentric to cartesian
%
%   r = GEOCENTOCART(t, geocen)
%   r = GEOCENTOCART(t, geocen3)
%
%   Input:
%     t the trixial ellipsoid object
%     geocen an n x 2 array of the geocentric coordinates [phic, lamc]
%     geocen3 an n x 3 array of the geocentric coordinates [phic, lamc, h]
%   Output:
%     r an n x 3 array of cartesian points
%
%   phic and lamc are measured in degrees.  With geocen (an n x 2 array), h is
%   assumed to be 0 so that r lies on the ellipsoid.
%
%   See also CARTTOGEOCEN, CART2TOCART

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = [cosd(geocen(:,1)) .* [cosd(geocen(:,2)), sind(geocen(:,2))], ...
       sind(geocen(:, 1))];
  r = r ./ vecabs(r./t.axes);
  if size(geocen, 2) == 3
    r = cart2tocart(t, r, geocen(:,3));
  end
end
