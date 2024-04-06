function [z, iz] = remx(x, y, alt)
%REMX   The remainder function
%
%   [z, iz] = REMX(x, y, alt)
%
% z is the remainder of x on division by y.  alt is optional flag, default
% false.  Result is in
%     [-y/2, y/2) if alt = false (the default)
%     (-y/2, y/2] if alt = true
% x can be any shape.  y should be a positive scalar.  iz corresponding
% integer quotient, s.t., x = iz * y + z.

  z = rem(x, y);
  z(z < -y/2) = z(z < -y/2) + y;
  z(z >= y/2) = z(z >= y/2) - y;
  if nargin == 3 && alt
    z(z == -y/2) = y/2;
  end
  iz = round((x - z) / y);
end
