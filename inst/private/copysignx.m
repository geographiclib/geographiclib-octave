function z = copysignx(x, y)
%COPYSIGNX   Copy the sign
%
%   COPYSIGNX(x,y) returns the magnitude of x with the sign of y.  x and y
%   can be any compatible shapes.
  z = abs(x) .* (1 - 2 * signbitx(y));
end
