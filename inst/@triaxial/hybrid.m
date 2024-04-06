function [pos2, dir2, s12, m12, M12, M21] = ...
      hybrid(t, pos1, dir1, cond, omgp, r2)
%HYBRID  the hybrid geodesic problem for a triaxial ellipsoid
%
%   [ellip2, alp2, s12, m12, M12, M21] = HYBRID(t, ellip1, alp1, cond)
%   [r2, v2, s12, m12, M12, M21] = HYBRID(t, r1, v1, cond)
%   [ellip2, alp2, s12, m12, M12, M21] = HYBRID(t, ellip1, alp1, cond, omgp,r2)
%   [r2, v2, s12, m12, M12, M21] = HYBRID(t, r1, v1, cond, omgp,r2)
%
%   Input:
%     t the trixial ellipsoid object
%     ellip1 1 x 2 vector of ellipsoidal starting coordinates [bet, omg]
%     alp1 the initial azimuths at ellip1
%     r1 1 x 3 vector for cartesian the starting point
%     v1 1 x 3 array of the starting direction at r1
%     cond k x 3 set of exit conditions
%     omgp (default 0), if set use [0, 360) as range for omg
%     r2 (default nan(1,3)), needed for cond(:,1) == 14
%   Output:
%     ellip2 1 x 2 vector of ellipsoidal coordinates at the endpoint
%     alp2 the azimuths at ellip2
%     r2 1 x 3 vector for cartesian endpoint
%     v2 1 x 3 vector for direction at r2
%     s12 the distance from ellip1 to ellip2 or r1 to r2
%     m12, M12, M21 the reduced length and geodesic scales
%
%   Solve the hybrid geodesic problem, namely given a point ellip1 on the
%   ellipsoid and an azimuth alp1, find the point ellip2 along the geodesic
%   satisfying one of the conditions cond.
%
%   cond is kx3 vector specifying k exit conditions
%     cond(:,1) = quantity tested
%        0 = s12 (but you should use cartdirect instead)
%        [1-10] any of the components of the ode state [r, v, m, mp, M, Mp]
%        [11-13] bet, omg, alp -- see also omgp flag
%        14 the directed distance to r2
%     cond(:,2) = value of corresponding variable
%     cond(:,3) = direction
%   Typical examples:
%     cond = [ 0, s12,  1] the direct problem
%     cond = [ 7, 0,   -1] find conjugate point
%     cond = [ 2, 0,   -1] find opposite umbilic point
%     cond = [ 2, 0,   -1] tracing a geodesic back to an umbilic point
%     cond = [11, bet2, 1;
%             13, 90,   1; % or [13, -90, -1]
%              7, 0,   -1] hybrid problem (2nd+3rd conds for bet2 = -bet1)
%     cond = [14, 0,   -1] closest approach to r2
%
%   The geodesic followed a maximum distance 2 * pi * t.a.  If the exit
%   condition has not been met, nans are returned.
%
%   This routine uses the function t.odesolver to solve the ODEs.  By default
%   this is @ode45 for Octave and @ode89 for MATLAB.  The default error
%   tolerances are multiplied by t.odemult which is set to 1e-10 by default.
%
%   Because this method is mainly for internal use by DISTANCE, this routine
%   is NOT vectorized.  Only a single starting point and directiom can be
%   specified.
%
%   See also CART2NORM, DISTANCE, RECKON, GEODDISTANCE, GEODRECKON

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  if nargin < 5, omgp = 0; end
  if nargin < 6, r2 = nan(1,3); end
  ellip = size(pos1, 2) == 2;
  if ellip
    [r1, v1] = elliptocart2(t, pos1, dir1);
  else
    [r1, v1] = deal(pos1, dir1);
  end
  l = cond(:, 1) <= 3 | cond(:, 1) == 6;
  cond(l, 2) = cond(l, 2) / t.b;
  [r2, v2, s12, m12, M12, M21] = ...
    hybridint(scaled(t), r1 / t.b, v1, cond, omgp, r2 / t.b);
  r2 = r2 * t.b; s12 = s12 * t.b; m12 = m12 * t.b;
  if ellip
    [pos2, dir2] = cart2toellip(t, r2, v2);
  else
    [pos2, dir2] = deal(r2, v2);
  end
end
