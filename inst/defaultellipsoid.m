function ellipsoid = defaultellipsoid(newdefault)
%DEFAULTELLIPSOID  Return or set the default ellipsoid
%
%   ellipsoid = DEFAULTELLIPSOID
%   DEFAULTELLIPSOID(newdefault)
%   DEFAULTELLIPSOID([])
%
%   The first form returns a vector of the equatorial radius and eccentricity
%   for the "default ellipsoid".  Initially the WGS84 ellispoid is the
%   default.  The second form allows a new default ellipsoid to be specified.
%   The third form restores the WGS84 ellipsoid as the default.  Use ecc2flat
%   and flat2ecc to convert between the eccentricity and the flattening.
%
%   See also ECC2FLAT, FLAT2ECC.

  narginchk(0, 1)
  persistent defellipsoid
  if nargin == 1
    assert(numel(newdefault) == 2 || numel(newdefault) == 0)
    defellipsoid = newdefault;
  end
  if isempty(defellipsoid)
    a = 6378137;
    f = 1/298.257223563;
    ecc = flat2ecc(f);
    defellipsoid = [a, ecc];
  end
  if nargin == 0
    ellipsoid = defellipsoid;
  end
end
