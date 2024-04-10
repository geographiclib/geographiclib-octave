function [ellip, alp, gam, latrad] = cart2toellip(t, r, v)
%CART2TOELLIP  Convert a surface point from cartesion to ellipsoidal
%
%   ellip = CART2TOELLIP(t, r)
%   [ellip, alp] = CART2TOELLIP(t, r, v)
%   [ellip, alp, gam, latrad] = CART2TOELLIP(t, r, v)
%
%   Input:
%     t the triaxial ellipsoid object
%     r an n x 3 array of cartesian points on the ellipsoid
%     v an n x 3 array of cartesian directions on the ellipsoid
%   Output:
%     ellip an n x 2 array of ellipsoidal coordinates [bet, omg]
%     alp an n x 1 array of azimuths alpha (only if v is given)
%     gam an n x 1 array of gamma, the invariant for a geodesic (if v is given)
%     latrad an n x 1 array of radii of curvature for a line of constant bet
%
%   bet, omg, and alp are measured in degrees.  This routine assumes that r lie
%   on the surface of the ellipsoid and that v is a unit vector tangent to the
%   ellipsoid at r.  To ensure that this is the case, call CARTNORM.  To
%   convert arbitrary points use CARTTOELLIP.
%
%   This routine can be called with r being a 1 x 3 array and v being an n x 3
%   array.  In this case, ellip is a 1 x 2 array and alp is an n x 1 array.
%
%   See also CARTNORM, CARTTOELLIP, ELLIPTOCART2, ELLIPTOCART

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  rn = r ./ t.axes;
  xi = rn(:,1); eta = rn(:,2); zeta = rn(:,3);
  g = vecdot([t.k2, t.k2 - t.kp2, -t.kp2], rn.^2);
  % Force umbilical point to be returned regardless of rounding errors in
  % multiplying and dividing by t.axes.
  g(vecabs(abs(r) - t.axes .* sqrt([t.kp2, 0, t.k2])) == 0) = 0;
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

  ellip = [bet, omg];

  if nargin < 2 || nargout < 2, return, end

  N = vecunit([-t.a * t.k2 * cb .* sb .* co ./ tx, ...
               -t.b * sb .* so, ...
                t.c * cb .* tz]);
  E = vecunit([-t.a * tx .* so, ...
                t.b * cb .* co, ...
                t.c * t.kp2 * sb .* co .* so ./ tz]);
  % At an oblate pole tx -> cb
  l = tx == 0;
  N(l,:) = [-co(l), -so(l), tx(l)] .* sb(l);
  E(l,:) = [-so(l),  co(l), tx(l)];
  % At a prolate pole tz -> so
  l = tz == 0;
  N(l,:) = [tz(l), -sb(l), cb(l)];
  E(l,:) = [tz(l),  cb(l), sb(l)] .* co(l);

  vn = vecdot(v, N);
  ve = vecdot(v, E);
  alp = atan2d(ve, vn);
  h = hypot(ve, vn);
  sa = ve./h; ca = vn./h;
  umb = cb == 0 & so == 0;
  sa(umb) = 1; ca(umb) = 0;

  % Calculate gamma with computed sines and cosines.
  gam = t.k2 * cb.^2 .* sa.^2 - t.kp2 * so.^2 .* ca.^2;
  tt = t.kp2 * so.^2;
  latrad = t.b * sqrt( (1 + t.e2 * tt) .* ...
                       (t.k2 * cb.^2 + tt) ./ ...
                       (t.k2 + tt) );

  % recompute cos(bet) and sin(omg) to see if r effectively corresponds to
  % an umbilical point.
  umb = cosd(bet) == 0 & sind(omg) == 0;
  if any(umb) && t.k2 ~= 0 && t.kp2 ~= 0
    w = sb .* co;
    up = vecunit(rn ./ t.axes);
    % return result in [-pi/2, pi/2]
    % 2*alp = atan2(s2, c2)
    % tan(2*alp) = 2*tan(alp)/(1-tan(alp)^2)
    % cos(alp) = sqrt((1+cos(2*alp))/2) = sin(2*alp) / sqrt(2*(1-cos(2*alp)))
    % sin(alp) = sqrt((1-cos(2*alp))/2) = sin(2*alp) / sqrt(2*(1+cos(2*alp)))
    alpu = atan2d(-v(:,2)./w / t.b, ...
                  (up(:,3).*v(:,1) - up(:,1).*v(:,3))./w * ...
                  (t.a * t.c / t.b^2)) / 2;
    % for northern umbilic flip by 180
    alpu(sb > 0) = alpu(sb > 0) - signx(alpu(sb > 0)) * 180;
    alp(umb) = alpu(umb);
  end
end
