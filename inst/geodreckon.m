function [lat2, lon2, azi2, S12, m12, M12, M21, a12_s12] = geodreckon ...
      (lat1, lon1, s12_a12, azi1, ellipsoid, flags)
%GEODRECKON  Point at specified azimuth, range on an ellipsoid
%
%   [lat2, lon2, azi2] = GEODRECKON(lat1, lon1, s12, azi1)
%   [lat2, lon2, azi2, S12, m12, M12, M21, a12_s12] =
%     GEODRECKON(lat1, lon1, s12_a12, azi1, ellipsoid, flags)
%
%   solves the direct geodesic problem of finding the final point and
%   azimuth given lat1, lon1, s12, and azi1.  The input arguments lat1,
%   lon1, s12, azi1, can be scalars or arrays of equal size.  lat1, lon1,
%   azi1 are given in degrees and s12 in meters.  The ellipsoid vector is
%   of the form [a, e], where a is the equatorial radius in meters, e is
%   the eccentricity.  If ellipsoid is omitted, the WGS84 ellipsoid (more
%   precisely, the value returned by defaultellipsoid) is used.  lat2,
%   lon2, and azi2 give the position and forward azimuths at the end point
%   in degrees.  The other outputs, S12, m12, M12, M21, a12 are documented
%   in geoddoc.  geoddoc also gives the restrictions on the allowed ranges
%   of the arguments.
%
%   flags (default 0) is a combination of 2 flags:
%      arcmode = bitand(flags, 1)
%      long_unroll = bitand(flags, 2)
%
%   If arcmode is unset (the default), then, in the long form of the call,
%   the input argument s12_a12 is the distance s12 (in meters) and the
%   final output variable a12_s12 is the arc length on the auxiliary sphere
%   a12 (in degrees).  If arcmode is set, then the roles of s12_a12 and
%   a12_s12 are reversed; s12_a12 is interpreted as the arc length on the
%   auxiliary sphere a12 (in degrees) and the corresponding distance s12 is
%   returned in the final output variable a12_s12 (in meters).
%
%   If long_unroll is unset (the default), then the value lon2 is in the
%   range [-180,180].  If long_unroll is set, the longitude is "unrolled"
%   so that the quantity lon2 - lon1 indicates how many times and in what
%   sense the geodesic encircles the ellipsoid.
%
%   The two optional arguments, ellipsoid and flags, may be given in any
%   order and either or both may be omitted.
%
%   When given a combination of scalar and array inputs, GEODRECKON behaves
%   as though the inputs were expanded to match the size of the arrays.
%   However, the setting up of the geodesic lines depends only on lat1 and
%   azi1.  So the combination of these arrays should be made as small as
%   possible.  For example, lat1, lon1, and azi1 could be specified as
%   column vectors of length M, while s12 is a row vector of length N.  The
%   result will be an M x N matrices however some portion of the compution
%   will performed on the smaller, size M, problem involving lat1 and azi1.
%   For N large, this speeds up the computation by about a factor of 2.
%
%   This is an implementation of the algorithm given in
%
%     C. F. F. Karney, Algorithms for geodesics,
%     J. Geodesy 87, 43-55 (2013);
%     https://doi.org/10.1007/s00190-012-0578-z
%     Addenda: https://geographiclib.sourceforge.io/misc/geod-addenda.html
%
%   This function duplicates some of the functionality of the RECKON
%   function in the MATLAB mapping toolbox.  Differences are
%
%     * When the ellipsoid argument is omitted, use the WGS84 ellipsoid.
%     * The routines work for prolate (as well as oblate) ellipsoids.
%     * The azimuth at the end point azi2 is returned.
%     * The solution is accurate to round off for abs(e) < 0.2.
%     * Redundant calculations are avoided when computing multiple
%       points on a single geodesic.
%     * Additional properties of the geodesic are calcuated.
%
%   See also GEODDOC, GEODDISTANCE, GEODAREA, DEFAULTELLIPSOID, FLAT2ECC.

% Copyright (c) Charles Karney (2012-2022) <charles@karney.com>.
%
% This is a straightforward transcription of the C++ implementation in
% GeographicLib and the C++ source should be consulted for additional
% documentation.  This is a vector implementation and the results returned
% with array arguments are identical to those obtained with multiple calls
% with scalar arguments.

  narginchk(4, 6)
  try
    % Use -0.0 for Z1 and Z to preserve the sign of 0.0.  geodreckon is one of
    % the few functions where this matters (to distinguish between
    % east-going and west-going meridional geodesics).
    S1 = size(lat1 + azi1); Z1 = -zeros(S1); % The size of geodesic line prob
    S = size(Z1 + lon1 + s12_a12); Z = -zeros(S); % The size of the full prob
  catch
    error('lat1, lon1, s12, azi1 have incompatible sizes')
  end
  if nargin <= 4
    ellipsoid = defaultellipsoid; flags = 0;
  elseif nargin == 5
    arg5 = ellipsoid(:);
    if length(arg5) == 2
      ellipsoid = arg5; flags = 0;
    else
      flags = arg5; ellipsoid = defaultellipsoid;
    end
  else
    arg5 = ellipsoid(:);
    arg6 = flags;
    if length(arg5) == 2
      ellipsoid = arg5; flags = arg6;
    else
      flags = arg5; ellipsoid = arg6;
    end
  end
  if length(ellipsoid) ~= 2
    error('ellipsoid must be a vector of size 2')
  end
  if ~isscalar(flags)
    error('flags must be a scalar')
  end

  arcmode = bitand(flags, 1);
  long_unroll = bitand(flags, 2);
  lat1 = AngRound(LatFix(lat1)) + Z1; azi1 = AngRound(azi1) + Z1;
  lat1 = lat1(:); azi1 = azi1(:);
  % expand lon1 to full size
  lon1 = lon1 + Z; lon1 = lon1(:);

  degree = pi/180;
  tiny = sqrt(realmin);

  a = ellipsoid(1);
  e2 = real(ellipsoid(2)^2);
  f = e2 / (1 + sqrt(1 - e2));
  f1 = 1 - f;
  ep2 = e2 / (1 - e2);
  n = f / (2 - f);
  b = a * f1;

  areap = nargout >= 4;
  redlp = nargout >= 5;
  scalp = nargout >= 6;

  A3x = A3coeff(n);
  C3x = C3coeff(n);

  [salp1, calp1] = sincosdx(azi1);
  [sbet1, cbet1] = sincosdx(lat1);
  sbet1 = f1 * sbet1; cbet1 = max(tiny, cbet1);
  [sbet1, cbet1] = norm2(sbet1, cbet1);
  dn1 = sqrt(1 + ep2 * sbet1.^2);

  salp0 = salp1 .* cbet1; calp0 = hypot(calp1, salp1 .* sbet1);
  ssig1 = sbet1; somg1 = salp0 .* sbet1;
  csig1 = cbet1 .* calp1; csig1(sbet1 == 0 & calp1 == 0) = 1; comg1 = csig1;
  [ssig1, csig1] = norm2(ssig1, csig1);

  k2 = calp0.^2 * ep2;
  epsi = k2 ./ (2 * (1 + sqrt(1 + k2)) + k2);
  A1m1 = A1m1f(epsi);
  C1a = C1f(epsi);
  B11 = SinCosSeries(true, ssig1, csig1, C1a);
  if ~arcmode
    s = sin(B11); c = cos(B11);
    stau1 = ssig1 .* c + csig1 .* s; ctau1 = csig1 .* c - ssig1 .* s;
    C1pa = C1pf(epsi);
  end
  C3a = C3f(epsi, C3x);
  A3c = -f * salp0 .* A3f(epsi, A3x);
  B31 = SinCosSeries(true, ssig1, csig1, C3a);

  if redlp || scalp
    A2m1 = A2m1f(epsi);
    C2a = C2f(epsi);
    B21 = SinCosSeries(true, ssig1, csig1, C2a);
  end

  if areap
    C4x = C4coeff(n);
    C4a = C4f(epsi, C4x);
    A4 = (a^2 * e2) * calp0 .* salp0;
    B41 = SinCosSeries(false, ssig1, csig1, C4a);
  end

  if prod(S) ~= prod(S1)
    if arcmode || abs(f) > 0.01
      C1a = expand(C1a, S1, Z);
    end
    k2 = expand(k2, S1, Z);
    dn1 = expand(dn1, S1, Z);
    A1m1 = expand(A1m1, S1, Z);
    B11 = expand(B11, S1, Z);
    C3a = expand(C3a, S1, Z);
    B31 = expand(B31, S1, Z);
    A3c = expand(A3c, S1, Z);

    salp0 = expand(salp0, S1, Z); calp0 = expand(calp0, S1, Z);
    ssig1 = expand(ssig1, S1, Z); csig1 = expand(csig1, S1, Z);
    somg1 = expand(somg1, S1, Z); comg1 = expand(comg1, S1, Z);
    if ~arcmode
      stau1 = expand(stau1, S1, Z); ctau1 = expand(ctau1, S1, Z);
      C1pa = expand(C1pa, S1, Z);
    end
    if redlp || scalp
      A2m1 = expand(A2m1, S1, Z);
      C2a = expand(C2a, S1, Z);
      B21 = expand(B21, S1, Z);
    end
    if areap
      % C4x isn't used again
      C4a = expand(C4a, S1, Z);
      A4 = expand(A4, S1, Z);
      B41 = expand(B41, S1, Z);
    end
  end

  if arcmode
    sig12 = s12_a12 * degree;
    [ssig12, csig12] = sincosdx(s12_a12);
    sig12 = sig12 + Z; sig12 = sig12(:);
    ssig12 = ssig12 + Z; ssig12 = ssig12(:);
    csig12 = csig12 + Z; csig12 = csig12(:);
  else
    s12_a12 = s12_a12 + Z; s12_a12 = s12_a12(:);
    tau12 = s12_a12 ./ (b * (1 + A1m1));
    s = sin(tau12); c = cos(tau12);
    B12 = - SinCosSeries(true, stau1 .* c + ctau1 .* s, ...
                         ctau1 .* c - stau1 .* s, C1pa);
    sig12 = tau12 - (B12 - B11);
    ssig12 = sin(sig12); csig12 = cos(sig12);
    if abs(f) > 0.01
      ssig2 = ssig1 .* csig12 + csig1 .* ssig12;
      csig2 = csig1 .* csig12 - ssig1 .* ssig12;
      B12 =  SinCosSeries(true, ssig2, csig2, C1a);
      serr = (1 + A1m1) .* (sig12 + (B12 - B11)) - s12_a12/b;
      sig12 = sig12 - serr ./ sqrt(1 + k2 .* ssig2.^2);
      ssig12 = sin(sig12); csig12 = cos(sig12);
    end
  end

  ssig2 = ssig1 .* csig12 + csig1 .* ssig12;
  csig2 = csig1 .* csig12 - ssig1 .* ssig12;
  dn2 = sqrt(1 + k2 .* ssig2.^2);
  if arcmode || redlp || scalp
    if arcmode || abs(f) > 0.01
      B12 = SinCosSeries(true, ssig2, csig2, C1a);
    end
    AB1 = (1 + A1m1) .* (B12 - B11);
  end
  sbet2 = calp0 .* ssig2;
  cbet2 = hypot(salp0, calp0 .* csig2);
  cbet2(cbet2 == 0) = tiny;
  somg2 = salp0 .* ssig2; comg2 = csig2;
  salp2 = salp0; calp2 = calp0 .* csig2;
  if long_unroll
    E = copysignx(1, salp0);
    omg12 = E .* (sig12 ...
                  - (atan2(   ssig2, csig2) - atan2(   ssig1, csig1)) ...
                  + (atan2(E.*somg2, comg2) - atan2(E.*somg1, comg1)));
  else
    omg12 = atan2(somg2 .* comg1 - comg2 .* somg1, ...
                  comg2 .* comg1 + somg2 .* somg1);
  end
  lam12 = omg12 + ...
          A3c .* ( sig12 + (SinCosSeries(true, ssig2, csig2, C3a) - B31));
  lon12 = lam12 / degree;
  if long_unroll
    lon2 = lon1 + lon12;
  else
    lon12 = AngNormalize(lon12);
    lon2 = AngNormalize(AngNormalize(lon1) + lon12);
  end
  lat2 = atan2dx(sbet2, f1 * cbet2);
  azi2 = atan2dx(salp2, calp2);
  if arcmode
    a12_s12 = b * ((1 + A1m1) .* sig12 + AB1);
  else
    a12_s12 = sig12 / degree;
  end
  a12_s12 = reshape(a12_s12, S);

  if redlp || scalp
    B22 = SinCosSeries(true, ssig2, csig2, C2a);
    AB2 = (1 + A2m1) .* (B22 - B21);
    J12 = (A1m1 - A2m1) .* sig12 + (AB1 - AB2);
    if redlp
      m12 = b * ((dn2 .* (csig1 .* ssig2) - dn1 .* (ssig1 .* csig2)) ...
                 - csig1 .* csig2 .* J12);
      m12 = reshape(m12, S);
    end
    if scalp
      t = k2 .* (ssig2 - ssig1) .* (ssig2 + ssig1) ./ (dn1 + dn2);
      M12 = csig12 + (t .* ssig2 - csig2 .* J12) .* ssig1 ./ dn1;
      M21 = csig12 - (t .* ssig1 - csig1 .* J12) .* ssig2 ./  dn2;
      M12 = reshape(M12, S); M21 = reshape(M21, S);
    end
  end

  if areap
    B42 = SinCosSeries(false, ssig2, csig2, C4a);
    salp12 = calp0 .* salp0 .* ...
             cvmgt(csig1 .* (1 - csig12) + ssig12 .* ssig1, ...
                   ssig12 .* ...
                   (csig1 .* ssig12 ./ max(1, (1 + csig12)) + ssig1), ...
                   csig12 <= 0);
    calp12 = salp0.^2 + calp0.^2 .* csig1 .* csig2;
    s = calp0 == 0 | salp0 == 0;
    salp12(s) = salp2(s) .* calp1(s) - calp2(s) .* salp1(s);
    calp12(s) = calp2(s) .* calp1(s) + salp2(s) .* salp1(s);
    if e2 ~= 0
      c2 = (a^2 + b^2 * eatanhe(1, e2) / e2) / 2;
    else
      c2 = a^2;
    end
    S12 = c2 * atan2(salp12, calp12) + A4 .* (B42 - B41);
    S12 = reshape(S12, S);
  end

  lat2 = reshape(lat2, S);
  lon2 = reshape(lon2, S);
  azi2 = reshape(azi2, S);

end

function B = expand(A, in, out)
% Expand a matrix to go from the setup of geodesic lines to the
% calculation of positions along the lines.
%
% A is a column vector of size prod(in).  in is a row vector representing
% the shape of A before it was converted to a column vector. out is a zero
% matrix with a shape of the full problem.  B is a column vector of size
% prod(size(out)).  in and out must be compatible, i.e.,
%
%   size( zeros(in) + out ) == size( out )
%
  k = size(A, 2);
  if k == 1
    B = reshape(A, in) + out;
    B = B(:);
  else
    % k > 1: Allow A/B to be k such columnized matrices
    B = repmat( out(:), 1, k );
    for i = 1 : k
      B(:, i) = expand(A(:, i), in, out);
    end
  end
end
