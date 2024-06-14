function [r2, v2, s12, m12, M12, M21] = hybridint(t, r1, v1, cond, altp, r2)
%HYBRIDINT  the internal version of the hybrid geodesic problem
%
%   [r2, v2, s12, m12, M12, M21] = HYBRIDINT(t, r1, v1, cond, omgp, r2)
%
%   This is the internal routine called by HYBRID.  This internal routine
%   assumed that t has been scaled so that t.b == 1.  And it only accepts
%   and returned positions and directions as cartesian coordinates.
%
%   See also HYBRID, DISTANCE, RECKON

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  if nargin < 5, altp = 0; end
  if nargin < 6, r2 = nan(1,3); end
  opt = odeset('AbsTol', 1e-6 * t.odemult, ...
               'RelTol', max(100*eps, 1.0e-3 * t.odemult), ...
               'Events', @(s, y) fexit(s, y, t, cond, altp, r2));
  % components of y are r,y,z, vx, vy, vz, m, m', M, M'
  y1 = [r1(:); v1(:); 0; 1; 1; 0];
  ie = [];
  if isfinite(sum(y1))
    tin = [0, 2 * 3.1416 * t.a]';          % Slightly more than 2*pi*a
    [~, ~, te, ye, ie] = t.odesolver(@(s, y) deriv(s, y, t), tin, y1, opt);
  end
  if isempty(ie)
    te = nan; ye = nan(1, 10);
  end
  r2 = ye(1, 1:3);
  v2 = ye(1, 4:6);
  s12 = te(1);
  m12 = ye(1, 7);
  M12 = ye(1, 9);
  M21 = ye(1, 8);
end

function yp = deriv(~, y, t)
% y and yp are column vectors.
  rr = y(1:3)'; vv = y(4:6)';
  [vv, acc, K] = accel(t, rr, vv);
  yp = [vv, acc]';

  if length(y) == 6, return; end

  m = y(7); mp = y(8); M = y(9); Mp = y(10);
  yp = [yp; mp; -K*m; Mp; -K*M];
end

function [w, term, dir] = fexit(~, y, t, cond, altp, rt)
  n = size(cond, 1);
  w = zeros(n, 1);
  term = ones(n, 1);
  dir = cond(:, 3);
  maxc = max(cond(:, 1));
  boa = [0, 0, 0];
  r = y(1:3)'; v = y(4:6)';
  [r, v] = cart2norm(t, r, v);
  y(1:6) = [r'; v'];
  if maxc > 10
    rr = y(1:3)';
    if maxc >= 13
      vv = y(4:6)';
      [boa(1:2), boa(3)] = cart2toellip(t, rr, vv);
      % if boa(3) < -90, boa(3) = boa(3) + 360; end
    else
      boa(1:2) = cart2toellip(t, rr);
    end
    if bitand(altp, 1) && boa(2) < 0, boa(2) = boa(2) + 360; end
    if bitand(altp, 2) && boa(2) < 0 && boa(1) > 0
      boa(1) = 180-boa(1); boa(2) = -boa(2);
    elseif bitand(altp, 4) && boa(2) > 0 && boa(1) > 0
      boa(1) = 180-boa(1); boa(2) = -boa(2);
    end
  end
  for i = 1:n
    c = cond(i, 1);
    switch c
      case 0
        w(i) = t;
      case {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
        w(i) = y(c);
      case {11, 12, 13}
        w(i) = boa(c - 10);
      case 14
        % goes from pos to neg as r passes closest to rt
        if vecdot(rt./t.axes.^2, r./t.axes.^2) > 0
          w(i) = vecdot(v, rt-r);
        else
          w(i) = vecabs(rt-r);
        end
    end
  end
  w = w - cond(:, 2);
end
