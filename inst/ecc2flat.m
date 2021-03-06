function f = ecc2flat(ecc)
%ECC2FLAT   Convert eccentricity to flattening
%
%   f = ECC2FLAT(ecc)
%
%   returns the flattening of an ellipsoid given the eccentricity.
%
%   See also FLAT2ECC.

  narginchk(1, 1)
  e2 = real(ecc.^2);
  f = e2 ./ (1 + sqrt(1 - e2));
end
