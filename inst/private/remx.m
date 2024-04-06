function z = remx(x, y)
%REMX   The remainder function
%
%   REMX(x, y) is the remainder of x on division by y.  Result is in [-y/2,
%   y/2].  x can be any shape.  y should be a positive scalar.

  z = rem(x, y);
  z(z < -y/2) = z(z < -y/2) + y;
  z(z >  y/2) = z(z >  y/2) - y;
end
