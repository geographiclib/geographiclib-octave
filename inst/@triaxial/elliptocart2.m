function [r, v] = elliptocart2(t, ellip, alp)
%ELLIPTOCART2  Convert a surface point from to ellipsoidal to cartesian
%
%   r = ELLIPTOCART2(t, ellip)
%   [r, v] = ELLIPTOCART2(t, ellip, alp)
%
%   Input:
%     t the triaxial ellipsoid object
%     ellip an n x 2 array of ellipsoid coordinates [bet, omg]
%     alp an n x 1 array of directions alpha
%   Output:
%     r an n x 3 array of cartesian coordinates
%     v an n x 3 array of cartesian directions
%
%   bet, omg, and alp are measured in degrees.  This routine can be called with
%   ellip being a 1 x 2 array and alp being an n x 1 array.  In this case, r is
%   a 1 x 3 array and v is an n x 3 array.  On the other hand, if ellip is an n
%   x 2 array and alp is a scalar, then both r and v are n x 3 arrays.
%
%   See also CARTNORM, CARTTOELLIP, ELLIPTOCART2, ELLIPTOCART

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  bet = ellip(:,1);
  omg = ellip(:,2);
  cb = cosd(bet); sb = sind(bet); co = cosd(omg); so = sind(omg);
  tx = sqrt(t.k2 * cb.^2 + t.kp2);
  tz = sqrt(t.k2 + t.kp2 * so.^2);
  r = [ t.a * co .* tx, ...
        t.b * cb .* so, ...
        t.c * sb .* tz ];

  if nargin < 3 || nargout < 2, return; end

  sa = sind(alp); ca = cosd(alp);
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

  v = ca .* N + sa .* E;

  umb = cb == 0 & so == 0;
  if any(umb) && t.k2 ~= 0 && t.kp2 ~= 0
    % sin(2*alp) and cos(2*alp)
    sa2 = 2 * sa(umb) .* ca(umb);
    ca2 = (ca(umb) - sa(umb)) .* (ca(umb) + sa(umb));
    % sign on 2nd component is -sign(cos(bet)*sin(omg)).  negative sign
    % gives normal conventipon of alpha measured clockwise.
    v(umb,:) = [t.a/sqrt(t.kp2) * co(umb) .* ca2, ...
                -t.b/sqrt(t.k2*t.kp2) * co(umb).* sb(umb) .* sa2, ...
                -t.c/sqrt(t.k2) * sb(umb) .* ca2];
  end
  v = vecunit(v);                      % Is this needed?
end
