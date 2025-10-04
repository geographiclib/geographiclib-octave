function [r2, h, count, err] = carttocart2(t, r)
%CARTTOCART2  Find the closest point on an ellipsoid
%
%   [R2, h] = CARTTOCART2(t, R)
%   [R2, h, count, err] = CARTTOCART2(t, R)
%
%   Input:
%     t the triaxial ellipsoid object
%     R an n x 3 array of cartesian points
%   Output:
%     R2 an n x 3 array of the closest cartesian points lying on t
%     h an n x 1 array of the directed distances, heights, from R2 to R
%     count an n x 1 array of the number of Newton's iterations
%     err an n x 1 array of the error in the function
%
%   See also CART2TOCART, CART2TOGEOD, CARTTOGEOD

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

% The formulation is due to Ligas, Stud. Geophys. Geod. (2012).  The major
% difference is to use an improved starting guess for Newton's method which
% guarantees convergence.  Panou and Korakitis, J. Geod. (2022) use the
% bisection method to solve the same equation; but this requires many more
% iterations to converge compared to Newton's method.  Some other points:
%
% * the method used here is simpler that PK (2022), in that there's no need
%   to treat special cases separately (indeed the method could easily be
%   generalized to deal with ellipsoids in higher dimensions);
%
% * a couple of tricks are used to improve the numerical accuracy of Newton's
%   method: (a) the function is evaluated with extended precision to avoid
%   the loss of precision; (b) the independent variable is offset so that the
%   singularity in the function is at the origin;
%
% * this routine stops at finding the closest point and defer computing
%   latitudes and longitudes for that to point to other functions.
  n = size(r, 1);
  ztol = t.b * eps / 8;
  tol = eps^(2/3);
  r(abs(r) <= ztol) = 0 * r(abs(r) <= ztol); % Preserve sign of +/-0
  count = zeros(n, 1);
  s = r .* t.axes;
  % axes2 = axes.^2 - c^2
  axes2 = (t.axes - t.c) .* (t.axes + t.c);
  p = max([    abs(s(:,  3))          , ...
            vecabs(s(:,2:3)) - axes2(2), ...
            vecabs(s       ) - axes2(1) ], [], 2 );
  l = count == 0;                       % Still to do
  c = ~l;                               % Convergence criteria met
  funp = @(p) f(t.axes, axes2, r, p);   % Capture t.axes, axes2, and r
  od = zeros(n, 1);
  for i = 1:50
    if ~any(l), break; end
    count(l) = count(l) + 1;
    [fv, fp] = funp(p);
    % We're done if f(p) <= 0 on initial guess; this can happens when z = 0.
    % However, since Newton converges from below, any negative f(p)
    % indicates convergence.
    l = l & ~(fv <= tol^2);
    d = -fv./fp;
    p(l) = p(l) + d(l);
    l = l & ~c;
    % converged if fv <= 8*eps (after first iteration) or
    % d <= max(sqrt(eps), |p|) * tol and d <= od.
    % N.B. d and od are always positive.
    c = d <= max(t.b^2 * sqrt(eps), p) * tol & d < od;
    if i > 1
      c = c | fv <= 8 * eps;
    end
    od = d;
  end
  err = funp(p);
  err(p == 0 & r(:, 3) == 0& err < 0) = 0;
  r2 = t.axes.^2 .* r./(p + axes2);
  % Deal with case p == 0 (when r2(3) is indeterminate).
  l = p == 0;
  if any(l)
    if axes2(1) == 0                     % sphere
      r2(l, 1) = r(l, 1);
    end
    if axes2(2) == 0                     % sphere or prolate
      r2(l, 2) = r(l, 2);
    end
    r2(l, 3) = t.axes(3) * r(l, 3) .* ...
        sqrt(1 - (r2(l, 1)/t.axes(1)).^2 + (r2(l, 2)/t.axes(2)).^2);
  end
  h = (p - t.axes(3)^2) .* vecabs(r2 ./ t.axes.^2);
end

function [fv, fp] = f(axes, axes2, r, p)
% In the general case, f has double poles at p = -axes2.  For p in [0, inf],
% f(p) decreases from inf to -1 (i.e., f'(p) is negative), and f''(p) is
% positive.  At p = c*abs(z), f(p) is positive.  Starting Newton's method at
% this point is guaranteed to converge.
%
% f(p) = (a*x/(p + a^2-c^2))^2 + (b*y/(p + b^2-c^2))^2 + (c*z/p)^2 - 1
% f(p) >= f1(p) where
% f1(p) = (c*z/p)^2 - 1
% f1(p) = 0 => p = c * |z|;
% f(p) >= f2(p) where
% f2(p) = (b*y/(p + b^2-c^2))^2 + (c*z/(p + b^2-c^2))^2 - 1
%       = (hypot(b*y, c*z)/(p + b^2-c^2))^2 - 1
% f2(p) = 0 => p = hypot(b*y, c*z) - (b^2-c^2)
% f(p) >= f3(p) where
% f3(p) = (a*x/(p+a^2-c^2))^2 + (b*y/(p+a^2-c^2))^2 + (c*z/(p+a^2-c^2))^2 - 1
%       = (hypot(a*x, b*y, c*z)/(p + a^2-c^2))^2 - 1
% f3(p) = 0 => p = hypot(a*x, b*y, c*z) - (a^2-c^2)
% lower bound on root of f(p) =
% max( c * |z|,
%      hypot(b*y, c*z) - (b^2-c^2),
%      hypot(a*x, b*y, c*z) - (a^2-c^2) )
%
% f(p) <= f4(p) where
% f4(p) = (a*x/p)^2 + (b*y/p)^2 + (c*z/p)^2 - 1
%       = (hypot(a*x, b*y, c*z)/p)^2 - 1
% f4(p) = 0 => p = hypot(a*x, b*y, c*z)
% upper bound on root of f(p) is hypot(a*x, b*y, c*z)
%
% [x0, y0, z0] = [a^2 * x/(p + a^2-c^2), b^2 * y/(p + b^2-c^2), c^2 * z/p]
%
% z0 is indeterminate if z = 0 and p = 0. x0 and y0 may also be
% indeterminate for x = 0 (for sphere) or y = 0 (for prolate).  Set these
% to zero.  Then determine z0 as signx(z) * c * (1 - (x/a)^2 - (y/b)^2).
  fp = 0*p; fv = fp - 1; fcorr = fp; d = eps/2;
  for k = 1:3
    l = r(:,k) ~= 0;
    g = (axes(k) * r(l,k) ./ (p(l) + axes2(k))).^2;
    % Accumulate fv is two pieces: multiples of d and the remainder.  The
    % first sums with no roundoff error for fv in [-1, 1].
    ga = round(g/d) * d; gb = g - ga;
    fv(l) = fv(l) + ga; fcorr(l) = fcorr(l) + gb;
    % fv(l) = fv(l) + g;
    fp(l) = fp(l) - 2 * g ./ (p(l) + axes2(k));
  end
  fv = fv + fcorr;
end
