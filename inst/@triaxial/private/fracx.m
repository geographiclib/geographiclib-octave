function [z, iz] = fracx(x, y, alt)
%FRACX the fraction part of x
%
%   [z, iz] = FRACX(x, y, alt)
%
% z is the remainder of x on division by y.  alt is optional flag, default
% false.  Result is in
%     [0, y) if alt = false (the default)
%     (0, y] if alt = true
% x can be any shape.  y should be a positive scalar.  iz corresponding
% integer quotient, s.t., x = iz * y + z.

  z = rem(x, y) + 0;                    % switch -0 to 0
  z(z <  0) = z(z <  0) + y;
  z(z >= y) = z(z >= y) - y;
  if nargin == 3 && alt
    z(z == 0) = y;
  end
  iz = round((x - z) / y);
end
