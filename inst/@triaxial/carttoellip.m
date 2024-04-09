function [ellip3, count, err] = carttoellip(t, r)
%CARTTOELLIP  Convert a general point from cartesian to ellipsoidal
%
%   ellip3 = CARTTOELLIP(t, r)
%   [ellip3, count, err] = CARTTOELLIP(t, r)
%
%   Input:
%     t the triaxial ellipsoid object
%     r an n x 3 array of cartesian points
%   Output:
%     ellip3 an n x 3 array of the ellipsoidal coordinates [bet, omg, u]
%     count an n x 1 array of the number of Newton's iterations
%     err an n x 1 array of the error in the function
%
%   bet and omg are measured in degrees.  u is the minor semiaxis of the
%   confocal ellipsoid on which r lies.  Thus, if r lies on t, u = t.c
%
%   See also ELLIPTOCART, CART2TOELLIP, ELLIPTOCART2

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  axes2 = (t.axes - t.c) .* (t.axes + t.c);
  ztol = t.b * eps / 8;
  r(abs(r) <= ztol) = 0 * r(abs(r) <= ztol); % Preserve sign of +/-0
  r2 = r.^2;

  [q, count, err] = newt(axes2, r2, t.b);

  axes = sqrt(axes2 + q);

  % This duplicates the code in cart2toellip with t.axes replaced by axes.  The
  % values of t.k2 and t.kp2 do not depend on u.

  rn = r ./ axes;
  xi = rn(:,1); eta = rn(:,2); zeta = rn(:,3);
  g = vecdot([t.k2, t.k2 - t.kp2, -t.kp2], rn.^2);
  % Force umbilic point to be returned regardless of rounding errors in
  % multiplying and dividing by axes.
  g(vecabs(abs(r) - axes .* sqrt([t.kp2, 0, t.k2])) == 0) = 0;
  h = sqrt(g.^2 + (4*t.k2*t.kp2) * eta.^2);
  cb = sqrt( (h + g)/2 ) / sqrt(t.k2);
  so = eta ./ cb;
  l = g < 0;
  so(l) = signx(eta(l)) .* sqrt( (h(l) - g(l))/2 ) / sqrt(t.kp2);
  cb(l) = abs(eta(l) ./ so(l));
  l = h == 0;
  so(l) = 0; cb(l) = 0;
  tz = sqrt(t.k2 + t.kp2 * so.^2);
  sb = zeta ./ tz;
  % at a prolate pole omg = 0 or 180, bet is arbitrary; pick bet = -90
  sb(tz == 0) = -1;
  tx = sqrt(t.k2 * cb.^2 + t.kp2);
  co = xi ./ tx;
  % at an oblate pole bet = +/-90, omg is arbirary; pick omg = 0
  co(tx == 0) = 1;
  bet = atan2d(sb, cb);
  omg = atan2d(so, co);

  ellip3 = [bet, omg, axes(:,3)];
end

function [t, count, err] = newt(axes2, r2, b)
% Solve
%   sum(r2(:,i) / (axes2(i) + t)) = 1
% for t.
%
% t is found by solving a cubic polynomial.  Since this solution might
% suffer from roundoff errors, the solution is refined using Newton's
% method.  The Newton's problem is very similar to that used in CARTTOCART2
% and, as in that routine, Newton's method is guaranteed to converge.

  tol = eps^(2/3);
  quadp = true;
  if ~quadp
    tmin = max([    r2(:,  3)               , ...
                    sum(r2(:,2:3), 2) - axes2(2), ...
                    sum(r2       , 2) - axes2(1)], [], 2);
    tmax = sum(r2, 2);
  else
    % Not worth doing with cubic starting pt
    tmin = max([quad(axes2(:,[2,3]),  r2(:,[2,3])), ...
                quad(axes2(:,[1,3]), [r2(:,1)+r2(:,2),r2(:,3)]), ...
                quad(axes2(:,[1,2]), [r2(:,1),r2(:,2)+r2(:,3)])], [], 2);
    tmax = min([quad(axes2(:,[2,3]), [r2(:,1)+r2(:,2),r2(:,3)]), ...
                quad(axes2(:,[1,3]), [r2(:,1),r2(:,2)+r2(:,3)])], [], 2);
  end
  n = size(r2, 1);
  count = zeros(n, 1);
  funp = @(t) f(axes2, r2, t);
  t = tmin;
  count = count + 1;
  fv = funp(t);
  % We're done if f(t) <= 0 on initial guess; this can happens when z = 0.
  l = ~(fv <= tol^2);                   % Still to do
  t3 = cubic(axes2, r2);
  t3 = max(tmin, min(tmax, t3));
  count(l) = count(l) + 1;
  [fv, fp] = funp(t3);
  t(l) = t3(l);
  l = l & ~(abs(fv) <= tol^2);
  d = -fv./fp;
  t(l) = max(tmin(l), t(l) + d(l));
  c = ~l;                               % Convergence criteria met
  od = zeros(n, 1);
  for i = 1:50
    if ~any(l), break; end
    count(l) = count(l) + 1;
    [fv, fp] = funp(t);
    % Since Newton converges from below, any negative f(t) indicates
    % convergence.
    l = l & ~(fv <= tol^2);
    d = -fv./fp;
    t(l) = t(l) + d(l);
    l = l & ~c;
    % converged if fv <= 8*eps or
    % d <= max(eps^(3/2), |t|) * tol and d <= od.
    % N.B. d and od are always positive.
    c = d <= max(b^2 * eps^(3/2), t) * tol & d < od;
    c = c | fv <= 8 * eps;
    %    i,fv,fp,d,od
    od = d;
  end
  err = funp(t);
  err(t == 0 & r2(:, 3) == 0 & err < 0) = 0;
end

function t = quad(l2, r2)
% Solve r2(1)/(t+l2(1)) + r2(2)/(t+l2(2)) - 1 = 0
  h = (sum(l2, 2) - sum(r2, 2)) / 2;
  c = prod(l2, 2) - l2(:,1).*r2(:,2) - l2(:,2).*r2(:,1);
  d = sqrt(max(0, h.^2 - c));
  t = d + abs(h);
  l = h > 0;
  t(l) = - c(l) ./ t(l);
end

function t = cubic(axes2, r2)
  [ea2, eb2] = deal(axes2(1), axes2(2));
  c = - ea2*eb2 * r2(:,3);
  b = ea2*eb2 - eb2 * r2(:,1) - ea2 * r2(:,2) - ...
      (ea2 + eb2) * r2(:,3);
  a = (ea2 + eb2) - sum(r2, 2);
  l = b > 0;
  % If b positsive there a cancellation in p = (3*b - a^2) / 3, so transform to
  % a polynomial in 1/t.  The resulting coefficients are
  [a(l), b(l), c(l)] = deal(b(l)./c(l), a(l)./c(l), 1./c(l));

  % Solve t^3 + a*t^2 + b*t + c = 0
  % see https://dlmf.nist.gov/1.11#iii
  p = (3*b - a.^2) / 3;
  q = (2*a.^3 - 9*a.*b + 27*c) / 27;

  % now switch to https://dlmf.nist.gov/4.43
  % We have 3 real roots, so 4*p^3 + 27*q^2 <= 0
  A = sqrt(max(0, -4/3 * p));
  ang = atan2(q, sqrt(max(0, -(4/27*p.^3 + q.^2))))/3;
  t = A/2 .* (cos(ang) * sqrt(3) - sin(ang)) - a / 3;

  t(l) = 1./t(l);
end

function [fv, fp] = f(axes2, r2, t)
  fp = 0*t; fv = fp - 1; fcorr = fp; d = eps/2;
  for k = 1:3
    l = r2(:,k) ~= 0;
    g = r2(l,k) ./ (t(l) + axes2(k));
    % Accumulate fv is two pieces: multiples of d and the remainder.  The
    % first sums with no roundoff error for fv in [-1, 1].
    ga = round(g/d) * d; gb = g - ga;
    fv(l) = fv(l) + ga; fcorr(l) = fcorr(l) + gb;
    % fv(l) = fv(l) + g;
    fp(l) = fp(l) - g ./ (t(l) + axes2(k));
  end
  fv = fv + fcorr;
end
