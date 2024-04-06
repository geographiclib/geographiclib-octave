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
