function z = atan2dx(y, x)
%ATAN2DX  Compute 2 argument arctangent with result in degrees
%
%   z = ATAN2DX(y, x) compute atan2(y, x) with result in degrees in
%   (-180,180] and quadrant symmetries enforced.  x and y must have
%   compatible shapes.

  persistent octavep
  if isempty(octavep)
    octavep = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  end
  if ~octavep
    % MATLAB implements symmetries already
    z = atan2d(y, x);
  else
    q1 = abs(y) > abs(x); q2 = ones(size(q1));
    x = q2 .* x; y = q2 .* y;           % expand x, y if necessary
    [ x(q1), y(q1) ] = deal( y(q1), x(q1) ); % swap
    q2 = signbitx(x); x(q2) = -x(q2);
    q = 2 * q1 + q2;
    % Now x >= 0 and x >= abs(y), so z is in [-45, 45]
    z = atan2d(y, x);
    % t = q == 0;              z(t) =    0 + z(t);
    t = q == 1; z(t) = copysignx(180, y(t)) - z(t);
    t = q == 2; z(t) =            90        - z(t);
    t = q == 3; z(t) =           -90        + z(t);
  end
end
