function t0 = scaled(t)
%SCALED  Return a scaled ellipsoid
%
%   t0 = t. SCALED
%
%   Return a similar ellipsoid with the middle semiaxes scaled to 1.  This
%   is used for integrating the ODEs for the geodesic to ensure all the
%   variables in the equations are of order unity.

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  t0 = t;
  if t0.b == 1, return; end
  t0.axes = t.axes/t.b;
  [t0.a, t0.b, t0.c] = deal(t.a/t.b, 1, t.c/t.b);
end
