function [sinx, cosx] = sincosdx(x)
%SINCOSDX  Compute sine and cosine with argument in degrees
%
%   [sinx, cosx] = SINCOSDX(x) compute sine and cosine of x in degrees with
%   exact argument reduction and quadrant symmetries enforced.

  persistent octavep
  if isempty(octavep)
    octavep = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  end
  if ~octavep
    % MATLAB implements argument reduction and symmetries already
    sinx = sind(x); cosx = cosd(x);
  else
    r = rem(x, 360);
    q = round(r / 90);
    r = r - 90 * q;
    q = mod(q, 4);
    r = r * (pi/180);
    sinx = sin(r); cosx = cos(r);
    t = q == 1; [sinx(t), cosx(t)] = deal( cosx(t), -sinx(t));
    t = q == 2; [sinx(t), cosx(t)] = deal(-sinx(t), -cosx(t));
    t = q == 3; [sinx(t), cosx(t)] = deal(-cosx(t),  sinx(t));
    sinx(sinx == 0) = copysignx(sinx(sinx == 0), x(sinx == 0));
    cosx = 0 + cosx;
  end
end
