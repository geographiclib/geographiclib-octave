function y = signx(x)
% return +/-1 depending on the signbit

  persistent octavep
  if isempty(octavep)
    octavep = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  end
  if octavep
    y = 1 - 2*signbit(x);
  else
    y = sign(x);
    y(x == 0) = sign(1 ./ x(x == 0));
  end
end
