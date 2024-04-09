function out = convert(t, in, from, to)
%CONVERT  convert coordinate between different representations
%
%   out = CONVERT(t, in, from, to)
%
%   Input:
%     in the input coordinates as n x 2 or n x 3 arrays
%     from the input format as a character string
%     to the output format as a character string
%   Output
%     out the output coordinates as n x 2 or n x 2 arrays
%
%   For arbitrary points, from and to can be
%     'cartesian'   [x, y, z] as general points
%     'geodetic'    [phi, lam, h]
%     'ellipsoidal' [bet, omg, u]
%
%   For points restricted to the surface of the ellipsoid, from and to can be
%     'cartesian2'  [x, y, z] confined to the surface of the ellipsoid
%     'geodetic'    [phi, lam]
%     'parametric'  [phip, lamp]
%     'geocentric'  [phic, lamc]
%     'ellipsoidal' [bet, omg]
%
%   If from = 'cartesian2', the input points are assumed to lie on the
%   ellipsoid.  If necessary, use CART2NORM to enforce this constraint.  An
%   error is thrown if to = 'cartesian2' or 'parametric' or 'geocentric' and
%   yet the input data are general 3d points.
%
%   See also CART2NORM, CART2TOELLIP, CARTTOELLIP, ELLIPTOCART, CART2TOGEOD,
%     CARTTOGEOD, GEODTOCART, CART2TOPARAM, PARAMTOCART2, CART2TOGEOCEN,
%     GEOCENTOCART2

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  if strcmp(from, to)
    out = in;
    return;
  end
  switch from
    case 'cartesian'
      threed = true;
      r = in;
    case 'cartesian2'
      threed = false;
      r = in;
    case 'geodetic'
      threed = size(in, 2) == 3;
      r = geodtocart(t, in);
    case 'parametric'
      threed = false;
      r = paramtocart2(t, in);
    case 'geocentric'
      threed = false;
      r = geocentocart2(t, in);
    case 'ellipsoidal'
      threed = size(in, 2) == 3;
      r = elliptocart(t, in);
    otherwise
      error(['unknown input type ', from])
  end
  % At this point the position is represented by r
  switch to
    case 'cartesian'
      out = r;
    case 'cartesian2'
      if threed
        error(['convert: cannot convert 3d points to ', to]);
      else
        out = r;
      end
    case 'geodetic'
      if threed
        out = carttogeod(t, r);
      else
        out = cart2togeod(t, r);
      end
    case 'parametric'
      if threed
        error(['convert: cannot convert 3d points to ', to]);
      else
        out = cart2toparam(t, r);
      end
    case 'geocentric'
      if threed
        error(['convert: cannot convert 3d points to ', to]);
      else
        out = cart2togeocen(t, r);
      end
    case 'ellipsoidal'
      if threed
        out = carttoellip(t, r);
      else
        out = cart2toellip(t, r);
      end
    otherwise
      error(['unknown output type ', to])
  end
end
