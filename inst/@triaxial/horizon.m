function r = horizon(t, ang, viewpt)
%HORIZON  Return a set of points on the horizon of the ellipsoid
%
%   R = HORIZON(t, ang, viewpt)
%
%   Input:
%     t the triaxial ellipsoid object
%     ang a vector of length n of angles defining the positions on the horizon
%     viewpt the geodetic coordinates defining the viewing direction
%   Output:
%     R an n x 3 array of 3d points on the horizon
%
%   ang is measured in degrees.  Typically this should include the full set of
%   angles, e.g., [-180:2:180] The viewpt is defined by the geodetic
%   coordinates [phi, lam].  The horizon is defined by the condition that the
%   normal on the horizon is orthogonal to normal at viewpt.
%
%   See also CARTPROJ

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

% Here's the theory.  Let u be the normal at viewpt.  The normal to the
% ellipsoid at a point r is r/axes^2.  The points on the horizon satisfy
% (r/axes^2) . u = 0; but this is equivalent to r/axes . u/axes = 0.  Now
% r/axes are just points on a sphere, and we therefore find these horizon
% points with the viewpoint view/axes.  Finally multiply the resulting
% spherical horizon points by axes to get the ellipsoidal horizon points.
  ang = ang(:);
  [phi, lam] = deal(viewpt(1), viewpt(2));
  view = [cosd(lam)*cosd(phi), sind(lam)*cosd(phi), sind(phi)];
  view1 = view ./ t.axes;
  lam1 = atan2d(view1(2), view1(1));
  phi1 = atan2d(view1(3), hypot(view1(2), view1(1)));
  M = rotm([phi1,lam1]); east = M(1,:); north = M(2,:);
  r = (cosd(ang)*east + sind(ang)*north) .* t.axes;
end
