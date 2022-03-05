function [sinx, cosx] = sincosde(x, t)
%SINCOSDX  Compute sine and cosine with argument in degrees
%
%   [sinx, cosx] = SINCOSDX(x, t) compute sine and cosine of x + t in
%   degrees with exact argument reduction and quadrant symmetries enforced.
%   x should be in [-180, 180] and t is a small correction

  q = round(x / 90);
  r = x - 90 * q;
  q = mod(q, 4);
  r = AngRound(r + t) * (pi/180);
  sinx = sin(r); cosx = cos(r);
  t = q == 1; [sinx(t), cosx(t)] = deal( cosx(t), -sinx(t));
  t = q == 2; [sinx(t), cosx(t)] = deal(-sinx(t), -cosx(t));
  t = q == 3; [sinx(t), cosx(t)] = deal(-cosx(t),  sinx(t));
  sinx(sinx == 0) = copysignx(sinx(sinx == 0), x(sinx == 0));
  cosx = 0 + cosx;
end
