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
%   from and to can be
%     'cartesian'   [x, y, z] as general points
%     'cartesian2'  [x, y, z] confined to the surface of the ellipsoid
%     'geodetic'    [phi, lam] or [phi, lam, h]
%     'parametric'  [phip, lamp] or [phip, lamp, h]
%     'geocentric'  [phic, lamc] or [phic, lamc, h]
%     'ellipsoidal' [bet, omg] or [bet, omg, u]
%
%   If from = 'cartesian2', the input points are assumed to lie on the
%   ellipsoid.  If necessary, use CART2NORM to enforce this constrainst.  An
%   error is thrown if to = 'cartesian2' and yet the input data are general 3d
%   points.
%
%   See also CART2NORM, CARTTOCART2, CART2TOCART, CART2TOELLIP,
%     CARTTOELLIP, ELLIPTOCART, CART2TOGEOCEN, GEOCENTOCART, CART2TOGEOD,
%     GEODTOCART, CART2TOPARAM, PARAMTOCART

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  if strcmp(from, to)
    out = in;
    return;
  end
  threed = size(in, 2) == 3;
  if threed, h = in(:, 3); end
  switch from
    case 'cartesian'
      switch to
        case 'ellipsoidal'
          out = carttoellip(t, in);
          return;
        otherwise
          [r, h] = carttocart2(t, in);
      end
    case 'cartesian2'
      r = in;
      threed = false;
      h = nan;
    case 'geodetic'
      r = geodtocart(t, in(:, 1:2));
    case 'parametric'
      r = paramtocart(t, in(:, 1:2));
    case 'geocentric'
      r = geocentocart(t, in(:, 1:2));
    case 'ellipsoidal'
      r = elliptocart(t, in);
      switch to
        case 'cartesian'
          out = r;
          return;
        otherwise
          if threed
            [r, h] = carttocart2(t, r);
          end
      end
    otherwise
      error(['unknown input type ', from])
  end
  % At this point the position is represented by r or, if threed, [r, h]
  switch to
    case 'cartesian'
      if threed
        out = cart2tocart(t, r, h);
      else
        out = r;
      end
      return;
    case 'cartesian2'
      if threed
        error('convert: cannot convert 3d points to cartesian2');
      else
        out = r;
      end
      return;
    case 'geodetic'
      out = cart2togeod(t, r);
    case 'parametric'
      out = cart2toparam(t, r);
    case 'geocentric'
      out = cart2togeocen(t, r);
    case 'ellipsoidal'
      if threed
        out = carttoellip(t, cart2tocart(t, r, h));
      else
        out = cart2toellip(t, r);
      end
      return;
    otherwise
      error(['unknown output type ', to])
  end
  if threed
    out = [out, h];
  end
end
