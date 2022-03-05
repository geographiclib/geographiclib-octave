function [d, t] = AngDiff(x, y)
%ANGDIFF  Compute angle difference accurately
%
%   [d, t] = ANGDIFF(x, y) computes z = y - x, reduced to (-180,180].  d =
%   round(z) and t = z - round(z).  x and y can be any compatible shapes.

  [d, t] = sumx(remx(-x, 360), remx(y, 360));
  [d, t] = sumx(remx(d, 360), t);
  l = d == 0 | abs(d) == 180;
  if any(l)
    z = y -x;
    d(l) = copysignx(d(l), cvmgt(z(l), -t(l), t(l) == 0));
  end
end
