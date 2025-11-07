classdef triaxial
%TRIAXIAL  A class describing triaxial ellipsoids
%
%   t = TRIAXIAL
%   t = TRIAXIAL(axes)
%   t = TRIAXIAL(ellipsoid)
%
%   Input:
%     axes = [a, b, c], the semiaxes
%     ellipsoid = [b, e2, k2, kp2], the ellipsoid parameters
%   Output:
%     t the triaxial ellipsoid object
%
%   The semiaxes can be in any units but must be ordered so that
%
%     inf > a >= b >= c > 0
%
%   The ellipsoid parameters are
%
%     e2  = (a^2 - c^2) / b^2
%     k2  = (b^2 - c^2) / (a^2 - c^2)
%     kp2 = (a^2 - b^2) / (a^2 - c^2)
%
%   The input parameters k2 and kp2 are scaled so that k2 + kp2 = 1.
%
%   If axes is omitted, the parameters of the triaxial reference ellipsoid
%   for the Earth, rounded to the nearest meter, are used:
%
%     a = 6378172, b = 6378102, c = 6356752
%
%   see
%     Panou, Korakitis, Pantazis, J. Geod. Sci. 10, 69 (2020)
%     https://doi.org/10.1515/jogs-2020-0105
%     Hu, Shum, Bevis, J. Geod. 97, 29 (2023)
%     https://doi.org/10.1007/s00190-023-01717-1
%
%   The properties of t that can be changed are
%
%     linesimptol the tolerance for line simplification in CARTPROJ
%     odesolver the name of the ODE solver
%     odemult the multiplier for the error tolerances for the ODE solver
%
%   The coordinate conversions are documented in
%
%     C. F. F. Karney,
%     Jacobi's solution for geodesics on a triaxial ellipsoid (2025).
%     https://arxiv.org/abs/arxiv:2511.01621
%
%   This paper uses a different method for solving the direct geodesic
%   problem; however, the solution using ODEs is discussed.
%
%   For a summary of the class, type
%     help @triaxial/doc
%
%   See also DOC, DEMO, TESTS

% Copyright (c) Charles Karney (2024-2025) <karney@alum.mit.edu>.

  properties
    axes
    a
    b
    c
    e2
    k2
    kp2
    linesimptol
    odemult
    odesolver
  end
  methods
    function obj = triaxial(axes)
    % TRIAXIAL.TRIAXIAL  the triaxial constructor
    %
    %   t = TRIAXIAL
    %   t = TRIAXIAL(axes)
    %   t = TRIAXIAL(ellipsoid)
    %
    %   Input:
    %     axes = [a, b, c], the semiaxes
    %     ellipsoid = [b, e2, k2, kp2], the ellipsoid parameters
    %   Output:
    %     t the triaxial ellipsoid object
    %
    %   The semiaxes can be in any units but must be ordered so that
    %
    %     inf > a >= b >= c > 0
    %
    %   The initialization with [b, e2, k2, kp2] specifies the ellipsoid in
    %   terms of
    %
    %     e2 = (a^2 - c^2) / b
    %     k2 = (b^2 - c^2) / (a^2 - c^2)
    %     kp2 = (a^2 - b^2) / (a^2 - c^2)
    %
    %   so that
    %
    %     a = b * sqrt(1 + e2*kp2)
    %     c = b * sqrt(1 - e2*k2)
    %
    %   N.B., k2 and kp2 are normalized to make k2 + kp2 = 1.
    %
    %   If axes is omitted, the parameters of the triaxial reference
    %   ellipsoid for the Earth, rounded to the nearest meter, are used:
    %
    %     a = 6378172, b = 6378102, c = 6356752
    %
    %   see
    %
    %     Panou, Korakitis, Pantazis, J. Geod. Sci. 10, 69 (2020)
    %     https://doi.org/10.1515/jogs-2020-0105
    %     Hu, Shum, Bevis, J. Geod. 97, 29 (2023)
    %     https://doi.org/10.1007/s00190-023-01717-1

      if nargin < 1
        axes = [6378172, 6378102, 6356752];
        % The conventional geographical longitude of the x axis is -14.94 deg
      end
      axes = axes(:)';
      assert(numel(axes) == 3 || numel(axes) == 4);
      if numel(axes) == 3
        assert(isfinite(axes(1)) && ...
               axes(1) >= axes(2) && axes(2) >= axes(3) && axes(3) > 0);
        obj.axes = axes;
        obj.a = obj.axes(1);
        obj.b = obj.axes(2);
        obj.c = obj.axes(3);

        s = (obj.a - obj.c) * (obj.a + obj.c);
        obj.e2 = s / obj.b^2;
        if s == 0
          % The sphere is a nonuniform limit, we can pick any values in [0,1]
          % s.t. k2 + kp2 = 1.  Here we choose to treat the sphere as an
          % oblate ellipsoid.
          obj.kp2 = 0; obj.k2 = 1 - obj.kp2;
        else
          obj.kp2 = (obj.a - obj.b) * (obj.a + obj.b)/s;
          obj.k2  = (obj.b - obj.c) * (obj.b + obj.c)/s;
        end
      else                              % numel(axes) == 4
        ksum = axes(3) + axes(4);
        assert(ksum > 0 && axes(3) >= 0 && axes(4) >= 0);
        obj.b = axes(1);
        obj.e2 = axes(2);
        obj.k2 = axes(3)/ksum;
        obj.kp2 = axes(4)/ksum;
        assert(obj.e2 * obj.k2 < 1);
        obj.a = obj.b * sqrt(1 + obj.e2 * obj.kp2);
        obj.c = obj.b * sqrt(1 - obj.e2 * obj.k2);
        obj.axes = [obj.a, obj.b, obj.c];
        assert(isfinite(obj.axes(1)) && ...
               obj.axes(1) >= obj.axes(2) && ...
               obj.axes(2) >= obj.axes(3) && obj.axes(3) > 0);
      end
      obj.linesimptol = 5e-4;
      obj.odemult = 1e-10;
      if is_octave
        % obj.odesolver = @ode23;  painfully slow
        obj.odesolver = @ode45;
      else
        % obj.odesolver = @ode23;
        % obj.odesolver = @ode45;
        % obj.odesolver = @ode78;
        obj.odesolver = @ode89;         % fast and accurate
        % obj.odesolver = @ode113;
      end
    end
    [r2, h, count, err]              = carttocart2(t, r)
    r                                = cart2tocart(t, r2, h)

    [ellip, alp, gam, latrad]        = cart2toellip(t, r, v)
    [r, v]                           = elliptocart2(t, ellip, alp)
    [ellip3, count, err]             = carttoellip(t, r)
    r                                = elliptocart(t, ellip)

    geod                             = cart2togeod(t, r)
    geod3                            = carttogeod(t, r)
    r                                = geodtocart(t, geod)

    param                            = cart2toparam(t, r)
    r                                = paramtocart2(t, param)

    geocen                           = cart2togeocen(t, r)
    r                                = geocentocart2(t, geocen)

    out                              = convert(t, in, from, to)

    [r, v]                           = cart2rand(t, n)

    [pos2, dir2, m12, M12, M21]      = reckon(t, pos1, dir1, s12)
    [s12, dir1, dir2, m12, M12, M21, count] = distance(t, pos1, pos2)
    [pos2, dir2, s12, m12, M12, M21] = hybrid(t, pos1, dir1, cond, altp, r2)

    [r, v]                           = cart2norm(t, r, v)
    t0                               = scaled(t)
    [r2f, r2b]                       = cartproj(t, r, viewpt, varargin)
    r                                = horizon(t, ang, viewpt)
  end
  methods (Access = private)
    [v, acc, K]                      = accel(t, r, v)
    [r2, v2, s12, m12, M12, M21]     = hybridint(t, r1, v1, cond, altp, r2)
  end
  methods (Static = true)
    [ellipn, alpn, flip]             = ellipnorm(ellip, alp, alt)
    [ellipn, alpn]                   = ellipflip(ellip, alp)
                                       demo(n)
                                       doc
                                       tests
  end
end
