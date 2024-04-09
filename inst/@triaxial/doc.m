function doc
%TRIAXIAL.DOC information on the TRIAXIAL class
%
%   TRIAXIAL.DOC
%
%   prints summary documentation on the TRIAXIAL class
%
%   t = TRIXIAL(axes) -- set up ellipsoid in terms of the semiaxes.
%
%   In function names amd function parameters
%   * cart or r refers to arbitary points in space expressed as [x, y, z];
%   * cart2 or r2 refers to points on the surface of the ellipsoid;
%   * geod specifies a point in geodetic coordinates [phi, lam];
%   * geod3 specifies a point in geodetic coordinates [phi, lam, h];
%   * param specifices a point in parametric coordinates [phip, lamp];
%   * param3 specifices a point in parametric coordinates [phip, lamp, h];
%   * geocen specifices a point in geocentric coordinates [phic, lamc];
%   * geocen3 specifices a point in geocentric coordinates [phic, lamc, h];
%   * ellip specifies a point in ellipsoidal coordinates [het, omg],
%   * ellip3 specifies a point in ellipsoidal coordinates [het, omg, u],
%       where u is the minor semiaxis of the confocal ellipsoid.
%
%   If the third element, h = height, of geod3, param3, or geocen3 is omitted,
%   it is taken as 0.  If the third element, u = minor semiaxis of the
%   confocal ellipsoid, of ellip3 is omitted, it is taken as the minor
%   semiaxis of the ellipsoid.
%
%   On the surface of an ellpsoid, specify a direction as a unit velocity v
%   in cartesian coordinates, or as an azimuth, alp, clockwise from a line
%   of constant ellipsoidal longitude, omg.
%
%   All angles are given in degrees.  All distances are given in the units
%   of the semiaxes.
%
%   Conversion routines
%
%     [r2, h] = t.CARTTOCART2(r)
%     r = t.CART2TOCART(r2, h)
%
%     [ellip, alp] = t.CART2TOELLIP(r, v)
%     [r, v] = t.ELLIPTOCART2(ellip, alp)
%     ellip3 = t.CARTTOELLIP(r)
%     r = t.ELLIPTOCART(ellip)
%
%     geod = t.CART2TOGEOD(r)
%     geod3 = t.CARTTOGEOD(r)
%     r = t.GEODTOCART(geod3)
%
%     param = t.CART2TOPARAM(r)
%     param3 = t.CARTTOPARAM(r)
%     r = t.PARAMTOCART(param3)
%
%     geocen = t.CART2TOGEOCEN(r)
%     geocen3 = t.CARTTOGEOCEN(r)
%     r = t.GEOCENTOCART(geocen3)
%
%     out = t.CONVERT(in, from, to)
%
%   Random points on the ellipsoid (and directions)
%     [r, v] = t.CART2RAND(n)
%
%   Geodesic calcultions
%     [pos2, dir2] = t.RECKON(pos1, dir1, s12)
%     [s12, dir1, dir2] = t.DISTANCE(pos1, pos2)
%     [r2, v2, s12] = t.HYBRID(r1, v1, cond)

%   Other
%     CART2NORM -- force points to lie on the surface
%     TRIAXIAL.ELLIPNORM -- reduce ellipoidal coords to convential ranges
%     SCALED -- return a scaled ellipsoid
%     TRIAXIAL.DEMO -- run demonstrations
%     TRIAXIAL.DOC -- this documentation
%     TRIAXIAL.TESTS -- the test suite
%     CARTPROJ -- plot a geodesic
%     HORIZON -- return point on the horizon

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  help triaxial/doc
end
