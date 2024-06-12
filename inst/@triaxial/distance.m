function [s12, dir1, dir2, m12, M12, M21, count] = distance(t, pos1, pos2)
%DISTANCE  the inverse geodesic problem on a triaxial ellipsoid
%
%   [s12, alp1, alp2] = DISTANCE(t, ellip1, ellip2)
%   [s12, v1, v2] = DISTANCE(t, r1, r2)
%   [s12, alp1, alp2, m12, M12, M21, count] = DISTANCE(t, ellip1, ellip2)
%   [s12, v1, v2, m12, M12, M21, count] = DISTANCE(t, r1, r2)
%
%   Input:
%     t the triaxial ellipsoid object
%     ellip1, ellip2 n x 2 arrays of ellipsoidal coordinates [bet, omg]
%     r1, r2 n x 3 arrays of 3d cartesian coordinates
%   Output:
%     s12 an n x 1 array of distances between point1 and point2
%     alp1, alp2 n x 1 arrays of azimuths at ellip1 and ellip2
%     v1, v2 n x 3 arrays of cartesian directions at r1 and r2
%     m12, M12, M21 n x 1 arrays of the reduced length and geodesic scales
%     count n x 1 array of the number of iterations to find the solution
%
%   Solve the inverse geodesic problem, namely given two points ellip1, ellip2
%   on the ellipsoid, find the shortest path s12 between them and the forward
%   azimuths alp1, alp2 at each point.
%
%   bet, omg, alp are measured in degrees.  The positions at the two endpoints
%   can be specified as either ellipsoidal coordinates or cartesian
%   coordinates.  In the latter case, the points must lie on the ellipsoid; if
%   necessary use CART2NORM to ensure this condition.  Both point1 and point2
%   must use the same representation.  point1 and point2 must be compatible
%   arrays (i.e., either of them can consist of a single row.
%
%   This routine is lamely vectorized.  The distance calculation for the n
%   pairs of point1 and point2 are handled independently.
%
%   For random pairs of points on a triaxial model for the Earth, the mean
%   running time for this routine is 4 s for Octave and 0.1 s for MATLAB.
%   (Octave is about 40 times slower than MATLAB.  This is due to the
%   different implementations of the ODE solver.)
%
%   See also CART2NORM, RECKON, HYBRID, GEODDISTANCE, GEODRECKON

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  Z = 0 * (pos1 + pos2);
  n = size(Z, 1);
  pos1 = pos1 + Z;
  pos2 = pos2 + Z;
  ellip = size(Z, 2) == 2;
  if ellip
    [pos1, ~, flip1] = triaxial.ellipnorm(pos1);
    [pos2, ~, flip2] = triaxial.ellipnorm(pos2);
    r1 = elliptocart2(t, pos1);
    r2 = elliptocart2(t, pos2);
  else
    [r1, r2] = deal(pos1, pos2);
  end
  r10 = r1 / t.b; r20 = r2 / t.b; t0 = scaled(t);
  s12 = zeros(n, 1); m12 = s12; M12 = s12; M21 = s12; count = s12;
  v1 = zeros(n, 3); v2 = v1;
  for k = 1:n
    [s12(k), v1(k,:), v2(k,:), m12(k), M12(k), M21(k), count(k)] = ...
        distanceint(t0, r10(k,:), r20(k,:));
  end
  s12 = s12 * t.b; m12 = m12 * t.b;
  if ellip
    [~, dir1] = cart2toellip(t, r1, v1);
    [~, dir2] = cart2toellip(t, r2, v2);
    [~, dir1(flip1)] = triaxial.ellipflip(pos1(flip1,:), dir1(flip1));
    [~, dir2(flip2)] = triaxial.ellipflip(pos2(flip2,:), dir2(flip2));
  else
    [dir1, dir2] = deal(v1, v2);
  end
end

function [s12, v1, v2, m12, M12, M21, count] = distanceint(t, r1, r2)
  check = false;
  M = eye(3);
  v1 = r2 - r1; s12 = vecabs(v1); m12 = s12; M12 = 1; M21 = 1; count = 0;

  % Deal with coincident, nearly coincident points, and nans
  if ~(s12 > sqrt(eps))
    if s12 == 0
      ellip = cart2toellip(t, (r1+r2)/2);
      [~, v1] = elliptocart2(t, ellip, 0);
    end
    [~, v2] = cart2norm(t, r2, v1);
    [~, v1] = cart2norm(t, r1, v1);
    return;
  end

  % Amount to move points by to ensure special cases get treated
  slop = 16*eps;

  rr = [r1; r2];
  l = abs(rr) <= slop;
  rr(l) = 0;
  ellip = cart2toellip(t, rr);
  bet = ellip(:,1);
  if t.k2 > 0 && t.kp2 > 0
    % Move point close to umbilics to umblics.
    p0 = t.axes .* sqrt([t.kp2, 0, t.k2]);
    l = vecabs(abs(rr) - p0) <= slop;
    rr(l,:) = [signx(rr(l,1)), 0*l(l), signx(rr(l,3))] .* p0;
  end

  rad = zeros(2,1);
  if t.kp2 == 0
    % Oblate, make point 1 closest to oblate pole
    rad = hypot(rr(:,1), rr(:,2));
    swapp = abs(rr(2,3)) > abs(rr(1,3)) | ...
            (abs(rr(2,3)) == abs(rr(1,3)) & rad(2) < rad(1));
  elseif t.k2 == 0
    % Prolate, make point 1 closest to prolate pole
    rad = hypot(rr(:,2), rr(:,3));
    swapp = abs(rr(2,1)) > abs(rr(1,1)) | ...
            (abs(rr(2,1)) == abs(rr(1,1)) & rad(2) < rad(1));
  else
    % Find point with max |bet|, if both have equal |bet|, pick the one with
    % largest |X| (smallest |omg| or |pi-omg|).
    swapp = ((abs(bet(2)) > abs(bet(1))) | ...
             (abs(bet(2)) == abs(bet(1)) & abs(rr(2,1)) > abs(rr(1,1))));
  end
  if swapp, rr = rr([2,1],:); rad = rad([2,1]); end

  if t.kp2 == 0
    % oblate, rotate about Z axis
    if rad(1) > 0
      k = 1;                            % move point 1 to omg = 0
    elseif rad(2) > 0
      k = 2;                            % point 1 on pole, move point 2
    else
      k = 0;                            % both points at poles
    end
    if k
      vp = rr(k,:) / rad(k);
      M(1:2,1:2) = [vp(1), -vp(2); vp(2), vp(1)];
    end
  elseif t.k2 == 0
    % prolate, rotate about Y axis to put point 2 on bet = 0
    k = 2;
    if rad(k) > 0
      vp = rr(k,:) / rad(k);
      M(2:3,2:3) = [vp(2), -vp(3); vp(3), vp(2)];
    end
  end
  rr = rr * M;
  l = abs(rr) <= slop;
  rr(l) = 0;

  % Reflect about X=0 if X1 < 0 || (X1 == 0 && X2 < 0).
  % Reflect about Y=0 if Y1 < 0 || (Y1 == 0 && Y2 < 0).
  % Reflect about Z=0 if Z1 > 0.
  flip = [rr(1,1) < 0 | (rr(1,1) == 0 & rr(2,1) < 0), ...
          rr(1,2) < 0 | (rr(1,2) == 0 & rr(2,2) < 0), ...
          rr(1,3) > 0];                 % Logical version
  flip = 1 - 2*flip;                    % +/-1 version
  rr = rr .* flip + 0;                  % convert -0 to +0
  M = M * diag(flip);
  ellip = cart2toellip(t, rr);
  [bet, omg] = deal(ellip(:,1), ellip(:,2));
  omg(omg == -180) = 180;
  [r1, r2] = deal(rr(1,:), rr(2,:));
  [bet1, bet2] = deal(bet(1), bet(2));
  [omg1, omg2] = deal(omg(1), omg(2));
  assert(bet1+bet2 <= 64*eps && bet2-bet1 >= -64*eps);% About 1.5 nm on earth
  bet2 = max(bet1, min(-bet1, bet2)); % Deals with roundoff

  % At this point we have
  % -90 <= bet1 <= bet2 <= -bet1; 0 <= omg1 <= 90
  assert(-90 <= bet1 && bet1 <= bet2 && bet2 <= -bet1 && -bet1 <= 90 && ...
         0 <= omg1 && omg1 <= 90);
  % if abs(bet2) = bet1 or prolate, then omg1 <= abs(omg2) <= 180-omg1
  assert(~(abs(bet2) == bet1 || t.k2 == 0) || ...
         omg1 <= abs(omg2) && abs(omg2) <= 180-omg1+eps(180));
  % For oblate: omg1 = 0, omg2 >= 0
  assert(~(t.kp2 == 0) || ...
           omg1 == 0 && omg2 >= 0);
  % For prolate: bet2 = 0
  assert(~(t.k2 == 0) || ...
         bet2 == 0);
  done = false;                         % A marker

  if bet1 == -90 && omg1 == 0
    % point 1 is an umbilical point at (bet,omg) = (-90,0)
    assert(r2(2) >= 0);
    if vecabs(r1 + r2) == 0
      % point 2 is opposite umbilic (bet,omg) = (90,180)
      % this handles prolate and oblate cases too
      r0 = [0, 1, 0];
      v0 = [-sqrt(t.kp2), 0, sqrt(t.k2)]; % Already normalized
      cond = [2, 0, -1];                % Check crossing Y = 0 plane
      [~, v2, sa] = hybridint(t, r0, v0, cond);
      s12 = 2*sa;
      % m12 = 0 for any geodesics between opposite umbilical points
      m12 = 0;
      % M12 = M21 = -1 only for this symmetric case going through [0, 1, 0].
      M12 = -1; M21 = -1;
      v1 = [1,-1,1] .* v2;
      done = true;
    elseif r2(2) > 0
      % General second point 2 (not on middle ellipse) with point 1 = umbilic
      alp2 = atan2d(sqrt(t.kp2)*sind(omg2), sqrt(t.k2)*cosd(bet2));
      [~, v2] = elliptocart2(t, [bet2, omg2], alp2);
      cond = [2, 0, -1];                % Check crossing Y = 0 plane
      [~, v1, s12, m12, M12, M21] = hybridint(t, r2, -v2, cond);
      v1 = -v1;
      [M12, M21] = deal(M21, M12);
      done = true;
    end
  end
  if r1(2) == 0 && r2(2) == 0 && ~done
    % Both points on middle ellipse.  This excludes the case of two
    % opposite umbilical points (already treated above).
    E = signx(r1(3) * r2(1) - r1(1) * r2(3)); % Better to call this an SE flag.
    [~, v1] = cart2norm(t, r1, -[1, 0, 1] * E);
    cond = [14, 0, -1];
    [~, v2, s12, m12, M12, M21] = hybridint(t, r1, v1, cond, 0, r2);
    if check
      [~, ~, s12a] = hybridint(t, r1, -v1, cond, 0, r2); %#ok<*UNRCH>
      assert(~(s12 - s12a > 128*slop));
    end
    done = m12 >= 0 | bet1 ~= -90 | bet2 ~= 90;
    % Only cases still to do are with bet1 = -90, bet2 = 90
    assert(done || (bet1 == -90 && bet2 == 90));
  end
  if ~done
    [s12, v1, v2, m12, M12, M21, count] = distanceint2(t, r1, r2);
  end
  [~, v1] = cart2norm(t, r1, v1);
  [~, v2] = cart2norm(t, r2, v2);
  v1 = v1 * M'; v2 = v2 * M';
  if swapp, [v1, v2, M12, M21] = deal(-v2, -v1, M21, M12); end
end

function [s12, v1, v2, m12, M12, M21, count] = distanceint2(t, r1, r2)
% input:
%   t = ellipsoid
%   r1, r2 are end points
% output:
%   s12 is the distance
%   v1, v2 are directions
%   m12, M12, M21 is the reduced length and geodesic scale
%   count is the number of iterations of Newton's method
% This deals with the generic case with a specific constraints assumed
% for r1 and r2, with
%
% -90 <= bet1 <= bet2 <= -bet1; 0 <= omg1 <= 90
% if abs(bet2) = bet1, then omg1 <= abs(omg2) <= 180-omg1
% For oblate: omg1 = 0, omg2 >= 0
% For prolate: bet2 = 0, omg1 <= abs(omg2) <= 180-omg1
%
% and the points aren't very close
% and neither point is an umbilical point (or pole)

% If r1(2) = 0, then either bet1 = -90 or omg1 = 0.

% Treat case where r2(2) = 0 specially.  Because of ordering of omg1, omg2
% check m12 for geodesic running around principal ellipse; accept if m12 >=
% 0.  Otherwise alp1 lies in (-90, 90); this only happens if bet1 = -90 and
% bet2 = 90.

% Now r2(2) > 0.

% If bet1 = -90, then there are two umbilical geodesics U[12] and U[34].
% with omgu[34] = 0, omgu[12] = 180.  alp1 lies in in (-90, 90)

% If omg1 = 0, then there are two umbilical geodesics U[14] and U[23] and
% omgu[14] = 0, check omgu[23] = 180.  alp1 lies in (0, 180).

% If r1 doesn't lie on Y = 0 then there are distinct four umbilical
% geodesics through r1.  Numbering the umblic geodesics U1 thru U4 for alp
% in 1, 2, 3, 4th quadrant (reckoning clockwise from 0) and defining omgu1
% thru omgu4 as the longitudes where bet = bet2.  For bet2 = bet1, omgu1 =
% omgu4 = omg1.  In general, omgu[14] in (0, 180), omgu[23] in (-180,0).

% If omg2 = one of the omguM, then pick alp1 = alpuM.  Otherwise if omg2
% in omguM and omgu(M+1), then alp1 in (alpuM, alpu(M+1)).  These ranges
% can be labeled E, S, W, N.

% For S, put alp1 in [0, 360).

% For E, put omg in [0,360) for omg comparisions.
% For S, W, N, then put omg in [-180,180) for omg comparisions.

% If omg2 > 0, compute omgu[14] only -- S is disallowed.  Compute omgu1
% first if omg2 > 90.

% If omg2 < 0, compute omgu[23] only -- N is disallowed.  Compute omgu2
% first if omg2 < -90

% For bet2 = bet1, omgu1 = omgu4 = omg1.  Replace alp1u = 90, alp4u = -90.
% For omg2 > 0, we have E (because omg2 > omg1) and alp1 is in (alp1u,
% alp2u) (no need to compute umbilical geodesics).  For omg2 < 0, compute
% omgu2 first if omg2 < -90.

  fzerop = false;
  check = false;
  slop = 16*eps;
  slopa = 180/pi * slop;
  ellip = cart2toellip(t, [r1; r2]);
  [bet1, omg1, bet2, omg2] = deal(ellip(1,1), ellip(1,2), ...
                                  ellip(2,1), ellip(2,2));
  bet2 = max(bet1, min(-bet1, bet2));   % Deals with roundoff
  done = false;
  count = 0;

  % Neither point at umbilical point (excludes poles for prolate/oblate)
  assert(~(abs(bet1) == 90 && (omg1 == 0 || abs(omg1) == 180)) && ...
         ~(abs(bet2) == 90 && (omg2 == 0 || abs(omg2) == 180)));

  % If both points lie on beta = 0 (this includes the case of both points
  % having the same value of beta on a prolate ellipsoid), then the two
  % points are connected by a geodesic running along a line of latitude.
  % However, it's not a shortest geodesic unless m12 >= 0.
  if bet1 == 0
    % Because omg1 <= abs(omg2) <= 180-omg1, the direction is given by
    E = signx(omg2);
    alp1 = E*90;
    [~, v1] = elliptocart2(t, [bet1, omg1], alp1);
    cond = [14, 0, -1];
    [~, v2, s12, m12, M12, M21] = hybridint(t, r1, v1, cond, 0, r2);
    done = m12 >= -128*slop;
    if check
      cond = [14, 0, -1];
      [~, ~, s12a] = hybridint(t, r1, -v1, cond, 0, r2);
      assert(~(s12 - s12a > 128*slop));
    end
  end
  if done, return; end

  % lists for 4 quadrants of alp1t
  alpu = atan2d(sqrt(t.kp2)*sind(omg1), sqrt(t.k2)*cosd(bet1));
  alpu = [alpu, 180-alpu, alpu-180, -alpu]';
  umb = (t.axes .* sqrt([ t.kp2, 0, t.k2])) .* [-1,1,1; 1,1,1; -1,1,1; 1,1,1];
  omgu = nan(4,1);

  range = [];
  cond = [];
  omgp = 0;
  if t.kp2 == 0
    % Oblate.
    if bet2 ~= bet1
      range = [0, 180];
      cond = [11, bet2, 1; 13, 90, 1];
    else                                % bet2 == bet1
      range = [0, 180];
      cond = [11, bet2, 1];
    end
  elseif t.k2 == 0
    % Prolate
    if omg2 > 0
      % if bet1 = bet2 = 0 and omg1 and omg2 are both positive, geodesic
      % is along line of constant beta which has been handled already.
      assert(bet1 ~= 0);
      range = [-90, 90];
    else
      range = [90, 270];
    end
    cond = [11, bet2, 1];
  elseif bet1 == -90
    range = [-90, 90];
    if bet2 == 90
      cond = [2, 0, -1];                  % crossing Y = 0 plane
    else
      cond = [11, bet2, 1];
    end
  elseif bet2 == bet1 && omg2 > 0
    range = [90, alpu(2)];
    cond = [11, bet2, 1];
    omgp = 1;
  elseif omg1 == 0
    range = [0, 180];
    cond = [11, bet2, 1; 13, 90, 1];
  elseif omg2 == 0
    range = [alpu(3), alpu(4)];
    cond = [11, bet2, 1; 13, -90, -1];
  elseif abs(omg2) == 180
    range = [alpu(1), alpu(2)];
    cond = [11, bet2, 1; 13, 90, 1];
    omgp = 1;
  else
    % That's the end of the cases where we can avoid computing umbilics
    % Last test in condu is to stop geodesic at umbilical point (needed if
    % bet2 is close to 90).
    condu = [11, bet2, 1; 14, 0, -1];
    if bet2 == bet1
      omgu(1) = omg1; omgu(4) = omg1;
      alpu(1) = 90; alpu(4) = -90;
    end
    if omg2 > 0                         % compute omgu[14]
      omgu(2) = -179; omgu(3) = -1;
      if omg2 > 90
        order = [1, 4, 3];
      else
        order = [3, 4, 1];
      end
    else                                % compute omgu[23]
      omgu(4) = 1; omgu(1) = 179;
      if omg2 < -90
        order = [1, 2, 3];
      else
        order = [3, 2, 1];
      end
    end
    for k = order
      l = mod(k, 4) + 1;
      for q = [k, l]
        if isnan(omgu(q))
          alp1 = alpu(q);
          [~, v1] = elliptocart2(t, [bet1, omg1], alp1);
          r2t = hybridint(t, r1, v1, condu, 0, umb(q,:));
          ellip2t = cart2toellip(t, r2t);
          [bet2t, omgu(q)] = deal(ellip2t(1), ellip2t(2)); %#ok<ASGLU>
          if check
            assert(abs(bet2t - bet2) < 512 * slopa);
          end
          done = abs(omg2 - omgu(q)) <= slopa;
          if done, break; end
        end
      end
      if done, break; end
      omgs = [omgu(k); omg2; omgu(l)];
      if k == 1
        omgs(omgs < 0) = omgs(omgs < 0) + 360;
      end
      % Replace < by <= in case the addition of 360 causes equality.  This is
      % probably caught by the fussy test on omg2 == omgu(q) above.
      if omgs(1) <= omgs(2) && omgs(2) <= omgs(3)
        % Good ranking omgs(1) < omg2(2) < omgs(3)
        switch k
          case 1
            cond = [11, bet2, 1; 13, 90, 1];
          case {2, 4}
            cond = [11, bet2, 1];
          case 3
            cond = [11, bet2, 1; 13, -90, -1];
        end
        range = [alpu(k), alpu(l)];
        omgp = k == 1;
        break;
      end
    end
  end
  if ~done
    assert(~isempty(range) && ~isempty(cond));
    if range(2) < range(1), range(2) = range(2) + 360; end
    if omgp && omg2 < 0, omg2 = omg2 + 360; end

    cond = [cond; 7, 0, -1];
    fun = @(alp1t) domgf(alp1t, t, cond, r1, bet1, omg1, omg2, omgp);
    if fzerop
      % Maybe narrow range slightly, since fzero evaluates fun at the
      % end points.
      % range = range + [1,-1] * 3*eps*180/pi;
      % alp1 = fzero(fun, range);
    else
      [alp1, count] = newtonx(fun, range, 180/pi);
    end
  end
  cond = [14, 0, -1];
  [~, v1] = elliptocart2(t, [bet1, omg1], alp1);
  [~, v2, s12, m12, M12, M21] = hybridint(t, r1, v1, cond, 0, r2);
end

function [domg, domgp] = domgf(alp1t, t, cond, r1, bet1, omg1, omg2, omgp)
  if nargin < 8, omgp = 0; end
  [~, v1] = elliptocart2(t, [bet1, omg1], alp1t);
  [r2, v2, ~, m12] = hybridint(t, r1, v1, cond, omgp);
  if size(cond, 1) > 0 && cond(1,1) == 2
    % If Y = 0 is the stopping condition, make this so
    r2(2) = cond(1,2);
  end
  [ellip2t, alp2, ~, rad2] = cart2toellip(t, r2, v2); omg2t = ellip2t(2);
  if omgp && omg2t < 0, omg2t = omg2t + 360; end
  domg = omg2t - omg2;
  domgp = m12 / ( cosd(alp2) * rad2 );
end

function [x, count] = newtonx(f, x0, scale)
% [y, yp] = f(x) returns a function and its derivative
% Find x, s.t. f(x) = 0.
% f must be monotonic increasing with f(x0(1)) < 0 and f(x0(2)) > 0.
% This is not vectorized.
  if nargin < 3, scale = 1; end
  maxcount = 100;
  xm = x0(1); xp = x0(2);
  xa = (xm + xp) / 2;
  count = 0;
  while count < maxcount
    count = count+1;
    [y, yp] = f(xa);
    if ~(abs(y) > eps * scale)
      x = xa;
      % 'zeroval'
      break
    end
    if y < 0
      xma = max(xm, xa);
      xpa = xp;
    else
      xpa = min(xp, xa);
      xma = xm;
    end
    dx = - y / yp;
    xb = xa + dx;
    if xb <= xm || xb >= xp || yp <= 0 || ~isfinite(yp)
      % 'bisect'
      xb = (xma + xpa) / 2;
      if xb == xma || xb == xpa
        x = xb;
        % 'zerorange'
        break;
      end
      dx = xb - xa;
    else
      % 'newtonx'
    end
    % range = xpa-xma
    if ~(abs(dx) > eps^0.75 * scale)
      x = xb;
      % 'zerodiff'
      break
    end
    % abs(dx)
    % diff = [xma, xb, xpa] - [xm, xa, xp]
    xa = xb;
    xm = xma; xp = xpa;
  end
end
