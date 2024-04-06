function r = elliptocart(t, ellip)
%ELLIPTOCART  Convert a general point from ellipsoidal to cartesian
%
%   r = ELLIPTOCART(t, ellip)
%   r = ELLIPTOCART(t, ellip3)
%
%   Input:
%     t the trixial ellipsoid object
%     ellip an n x 2 array of the ellipsoidal coordinates [bet, omg]
%     ellip3 an n x 3 array of the ellipsoidal coordinates [bet, omg, u]
%   Output:
%     r an n x 3 array of cartesian points
%
%   bet and omg are measured in degrees.  With ellip (an n x 2 array), u is
%   assumed to be t.c so that r lies on the ellipsoid.  In this case,
%   ELLIPTOCART2 is called.
%
%   See also CARTTOELLIP, CART2TOELLIP, ELLIPTOCART2

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  if size(ellip, 2) ~= 3
    r = elliptocart2(t, ellip);
    return;
  end
  bet = ellip(:,1);
  omg = ellip(:,2);
  u = ellip(:,3);
  axes = sqrt(((t.axes - t.c) .* (t.axes + t.c)) + u.^2);
  axes(:, 3) = u;

  cb = cosd(bet); sb = sind(bet); co = cosd(omg); so = sind(omg);
  tx = sqrt(t.k2 * cb.^2 + t.kp2);
  tz = sqrt(t.k2 + t.kp2 * so.^2);
  r = [ co .* tx, cb .* so, sb .* tz ];
  r = axes .* r;
end
