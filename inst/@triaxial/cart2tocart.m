function r = cart2tocart(t, r2, h)
%CART2TOCART  Convert a surface point + height to a cartesian point
%
%   r = CART2TOCART(t, r2, h)
%
%   Input:
%     t the triaxial ellipsoid object
%     r2 an n x 3 array of cartesian points on the ellipsoid
%     h an n x 1 array of heights
%   Output:
%     r an n x 3 array of cartesian points a height h above r2
%
%   See also CARTTOCART2, CART2TOGEOD, CARTTOGEOD

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  r = r2 + h .* vecunit(r2 ./ t.axes.^2);
end
