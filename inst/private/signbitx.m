function z = signbitx(x)
%SIGNBITX   Copy the sign
%
%   SIGNBITX(x) returns true for the elements of x which have the signbit
%   set.  x can be any shape.

  persistent octavep
  if isempty(octavep)
    octavep = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  end
  if octavep
    z = signbit(x);
  else
    l = x < 0 | (x == 0 & 1./x < 0);
  end
end
