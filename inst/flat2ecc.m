function ecc = flat2ecc(f)
%FLAT2ECC   Convert flattening to eccentricity
%
%   ecc = FLAT2ECC(f)
%
%   returns the eccentricity of an ellipsoid given the flattening.
%
%   See also ECC2FLAT.

  narginchk(1, 1)
  ecc = sqrt(f .* (2 - f));
end
