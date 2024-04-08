function [pos2, dir2, m12, M12, M21] = reckon(t, pos1, dir1, s12)
%RECKON  the direct geodesic problem for a triaxial ellipsoid
%
%   [ellip2, alp2] = RECKON(t, ellip1, alp1, s12)
%   [r2, v2] = RECKON(t, r1, v1, s12)
%   [ellip2, alp2, m12, M12, M21] = RECKON(t, ellip1, alp1, s12)
%   [r2, v2, m12, M12, M21] = RECKON(t, r1, v1, s12)
%
%   Input:
%     t the trixial ellipsoid object
%     ellip1 n x 2 array of ellipsoidal coordinates [bet, omg]
%     alp1 n x 1 array of azimuths at ellip1
%     r1 n x 3 array of cartesian starting points
%     v1 n x 3 array of cartesian directions at r1
%     s12 an m long vector of distances
%   Output:
%     ellip2 m*n x 2 arrays of ellipsoidal coordinates [bet, omg]
%     alp2 m*n x 1 array of azimuths at ellip2
%     r2 an m*n x 3 array of cartesian positions at distances s12
%     v2 an m*n x 3 array of cartesian directions at r2
%     m12, M12, M21 m*n x 1 arrays of reduced length and geodesic scales
%
%   Solve the direct geodesic problem, namely given a point ellip1 on the
%   ellipsoid and an azimuth alp1, find the point ellip2 a distance s12 away
%   and the azimuth at that point.
%
%   Compute positions at distances s12 along the geodesic starting at point1.
%   bet, omg, alp are measured in degrees.  The starting point and direction
%   can be specified as either ellipsoidal coordinates or cartesian
%   coordinates.  In the latter case, r1 must lie on the ellipsoid and v1 must
%   be a unit vector tangent to the ellipsoid at r1; if necessary use
%   CART2NORM to ensure this condition.  The distances can be given in any
%   order.  Internally the positive and negative distances are handled by two
%   separated calls to the ODE solver.  You can "punctuate" the results by
%   inserting nans into s12.
%
%   This routine uses the function t.odesolver to solve the ODEs.  By default
%   this is @ode45 for Octave and @ode89 for MATLAB.  The default error
%   tolerances are multiplied by t.odemult which is set to 1e-10 by default.
%
%   This routine is NOT vectorized.  The positions along each geodesic are
%   calculated handled independently.  For each geodesic, the positions at all
%   the values of s12 are computed.  Thus ellip2(1:m, :) and alp2(1:m) are the
%   positions and azimuths at distances s12(1:m) along the geodesic starting
%   at ellip1(1,:) and alp1(1).
%
%   The cost of the computation is proportional to the maximum positive value
%   of s12 plus the maximum absolute negative value of s12.  There's no big
%   penalty in adding intermediate points along the geodesic.
%
%   If none of m12, M12, M21 are returned, a smaller set of differential
%   equations is solved.  This results in very slightly different results
%   returned for ellip2 and alp2.
%
%   See also CART2NORM, DISTANCE, HYBRID, GEODDISTANCE, GEODRECKON

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  ellip = size(pos1, 2) == 2;
  if ellip
    [r1, v1] = elliptocart2(t, pos1, dir1);
  else
    [r1, v1] = deal(pos1, dir1);
  end
  r1 = r1 / t.b; s12 = s12(:) / t.b;
  Z = 0 * (r1 + v1);
  n = size(Z, 1);
  r1 = r1 + Z;
  v1 = v1 + Z;
  m = length(s12);
  r2 = zeros(m*n, 3); v2 = r2;
  if nargout <= 2
    for k = 1:n
      ind = m * (k-1) + (1:m);
      [r2(ind,:), v2(ind, :)] = reckonint(scaled(t), r1(k,:), v1(k,:), s12);
    end
  else
    m12 = zeros(m*n, 1); M12 = m12; M21 = m12;
    for k = 1:n
      ind = m * (k-1) + (1:m);
      [r2(ind, :), v2(ind, :), m12(ind), M12(ind), M21(ind)] = ...
          reckonint(scaled(t), r1, v1, s12);
    end
    m12 = m12 * t.b;
  end
  r2 = r2 * t.b;
  if ellip
    [pos2, dir2] = cart2toellip(t, r2, v2);
  else
    [pos2, dir2] = deal(r2, v2);
  end
end

function [r2, v2, m12, M12, M21] = reckonint(t, r1, v1, s12)
  n = length(s12);
  if ~isfinite(sum([r1, v1])) || n == 0
    r2 = nan(n, 3); v2 = nan(n, 3);
    m12 = nan(n, 1); M12 = nan(n, 1); M21 = nan(n, 1);
    return
  end
  s12(~isfinite(s12)) = nan;            % Convert +/-inf to nan
  redlen = nargout > 2;                 % Compute reduced length
  y1 = [r1, v1]';
  if redlen
    vlen = 10;
    y1 = [y1; 0; 1; 1; 0];
  else
    vlen = 6;
  end
  mult = t.odemult;
  opt = odeset('AbsTol', 1e-6*mult, ...
               'RelTol', max(100*eps, 1.0e-3*mult));
  ode = t.odesolver;
  [s12, ~, ind] = unique(s12);
  s12m = s12(s12 < 0); nm = length(s12m); s12m = s12m(end:-1:1);
  nz = sum(s12 == 0);
  s12p = s12(s12 > 0); np = length(s12p);
  nnan = sum(isnan(s12));
  ntot = length(s12);
  assert(ntot == nm + nz + np + nnan);
  pts = nan(ntot, vlen);                % Preallocate result array
  if nz, pts(nm + 1, :) = y1'; end
  for i = 0:1
    k = 2;
    if i == 0
      if nm == 0, continue, end
      if nm == 1
        % This is a ode solver quirk; need at least three elements of time
        % variable to get same time points returned.
        tin = [0; s12m(1)/2; s12m];
        k = 3;
      else
        tin = [0; s12m];
      end
    else
      if np == 0, continue; end
      if np == 1
        % This is a ode solver quirk; need at least three elements of time
        % variable to get same time points returned.
        tin = [0; s12p(1)/2; s12p];
        k = 3;
      else
        tin = [0; s12p];
      end
    end
    [~, ptsa] = ode(@(s, y) deriv(s, y, t), tin, y1, opt);
    if i == 0
      pts(1:nm, :) = ptsa(end:-1:k,:);
    else
      pts(nm + nz + (1:np), :) = ptsa(k:end,:);
    end
  end
  ptsa = pts(ind,:);
  r2 = ptsa(:, 1:3);
  v2 = ptsa(:, 4:6);
  if redlen
    m12 = ptsa(:, 7); M12 = ptsa(:, 9);
    M21 = ptsa(:, 8);                   % AG Eq 29: dm12/ds2 = M21
  end

end

% Expression for acceleration (Panou 2019)
% -[r,y,z]/[a2,b2,c2]/(r2/a2^2+y2/b2^2+z2/c2^2)*(vx2/a2+vy2/b2+vz2/c2)
function yp = deriv(~, y, t)
% y and yp are column vectors.
  [v, acc, K] = accel(t, y(1:3)', y(4:6)');
  yp = [v, acc]';

  if length(y) == 6, return; end

  m = y(7); mp = y(8); M = y(9); Mp = y(10);
  yp = [yp; mp; -K*m; Mp; -K*M];
end
