function y = AngNormalize(x)
%ANGNORMALIZE  Reduce angle to range (-180, 180]
%
%   x = ANGNORMALIZE(x) reduces angles to the range (-180, 180].  x can be
%   any shape.

  y = remx(x, 360);
  y(abs(y) == 180) = copysignx(180, x(abs(y) == 180));
end
