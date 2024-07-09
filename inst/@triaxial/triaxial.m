classdef triaxial
%TRIAXIAL  A class describing triaxial ellipsoids
%
%   t = TRIAXIAL
%   t = TRIAXIAL(axes)
%
%   Input:
%     axes = [a, b, c] a length 3 vector of semiaxes
%   Output:
%     t the triaxial ellipsoid object
%
%   The semiaxes can be in any units but must be ordered so that
%
%     inf > a >= b >= c > 0
%
%   If axes is omitted, the parameters of the triaxial reference ellipsoid for
%   the Earth are used; see Hu, Shum, Bevis, J. Geod. 97, 29 (2023)
%   https://doi.org/10.1007/s00190-023-01717-1
%
%   The properties of t that can be changed are
%
%     linesimptol the tolerance for line simplification in CARTPROJ
%     odesolver the name of the ODE solver
%     odemult the multiplier for the error tolerances for the ODE solver
%
%   The case a == b == c. the sphere, is a nonuniform limit.  How
%   ellipsoidal coordinates in this case depends on the properties
%
%     k2 = (b^2 - c^2) / (a^2 - c^2)
%     kp2 = (a^2 - b^2) / (a^2 -  c^2)
%
%   By default, this class treats the sphere as a limiting form of an oblate
%   ellipsoid, i.e., k2 = 1, kp2 = 0.  However, these values can be
%   overridden.  The values must satisfy k2 >= 0, kp2 >= 0, k2 + kp2 == 1.
%   This choice governs the definition of ellipsoidal coordinates and
%   affects the inner workings of DISTANCE
%
%   See also DOC, DEMO, TESTS, CART2TOELLIP, DISTANCE

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

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
    %
    %   Input:
    %     axes a length 3 vector of semiaxes [a, b, c]
    %   Output:
    %     t the triaxial ellipsoid object
    %
    %   The axes can be in any units but must be ordered so that
    %
    %     inf > a >= b >= c > 0
    %
    %   If axes is omitted, the parameters of the triaxial reference
    %   ellipsoid for the Earth, rounded to the nearest meter, are used;
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
      assert(numel(axes) == 3 && isfinite(axes(1)) && ...
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
