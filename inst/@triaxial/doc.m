function doc
%TRIAXIAL.DOC information on the TRIAXIAL class
%
%   TRIAXIAL.DOC
%
%   prints summary documentation on the TRIAXIAL class
%
%   t = TRIAXIAL([a, b, c]) -- set up ellipsoid in terms of the semiaxes.
%   t = TRIAXIAL([b, e2, k2, kp2]) -- or in terms of ellipsoid parameters.
%
%   In function names and function parameters
%   * cart or R refers to arbitrary points in space expressed as [X, Y, Z];
%   * cart2 or R2 refers to points on the surface of the ellipsoid;
%   * geod specifies a point in geodetic coordinates [phi, lam];
%   * param specifies a point in parametric coordinates [phip, lamp];
%   * geocen specifies a point in geocentric coordinates [phic, lamc];
%   * ellip specifies a point in ellipsoidal coordinates [het, omg],
%   * geod3 specifies a point in geodetic coordinates [phi, lam, h];
%   * ellip3 specifies a point in ellipsoidal coordinates [het, omg, u],
%       where u is the minor semiaxis of the confocal ellipsoid.
%
%   If the third element, h = height, of geod3 is omitted, it is taken as 0.
%   If the third element, u = minor semiaxis of the confocal ellipsoid, of
%   ellip3 is omitted, it is taken as the minor semiaxis of the ellipsoid.
%
%   On the surface of an ellipsoid, specify a direction as a unit velocity V
%   in cartesian coordinates, or as an azimuth, alp, clockwise from a line
%   of constant ellipsoidal longitude, omg.
%
%   All angles are given in degrees.  All distances are given in the units
%   of the semiaxes.
%
%   The accuracy of the geodesic routines for can be found using the
%   test set available at
%
%     https://10.5281/zenodo.1251079
%
%   This consists of 500000 shortest geodesics on the ellipsoid uses
%   triaxial(sqrt([2, 1, 1/2])).  Scaling the errors by 6.4e6 m to obtain
%   the results for an ellipsoid roughly the size of the earth,
%   we have the following estimates of the maximum error:
%
%               MATLAB  Octave
%               R2023a  8.4.0
%
%     reckon    0.1 um  0.5 um
%     distance   10 um  500 m
%
%   Conversion routines
%
%     [R2, h] = t.CARTTOCART2(R)
%     R = t.CART2TOCART(R2, h)
%
%     [ellip, alp] = t.CART2TOELLIP(R, V)
%     [R, V] = t.ELLIPTOCART2(ellip, alp)
%     ellip3 = t.CARTTOELLIP(R)
%     R = t.ELLIPTOCART(ellip)
%
%     geod = t.CART2TOGEOD(R)
%     geod3 = t.CARTTOGEOD(R)
%     R = t.GEODTOCART(geod3)
%
%     param = t.CART2TOPARAM(R)
%     R = t.PARAMTOCART2(param)
%
%     geocen = t.CART2TOGEOCEN(R)
%     R = t.GEOCENTOCART2(geocen)
%
%     out = t.CONVERT(in, from, to)
%
%   Random points on the ellipsoid (and directions)
%     [R, V] = t.CART2RAND(n)
%
%   Geodesic calculations
%     [pos2, dir2] = t.RECKON(pos1, dir1, s12)
%     [s12, dir1, dir2] = t.DISTANCE(pos1, pos2)
%     [R2, V2, s12] = t.HYBRID(R1, V1, cond)
%
%   Other
%     CART2NORM -- force points to lie on the surface
%     TRIAXIAL.ELLIPNORM -- reduce ellipoidal coords to standard ranges
%     TRIAXIAL.ELLIPFLIP -- switch ellipoidal coords to the other sheet
%     SCALED -- return a scaled ellipsoid
%     TRIAXIAL.DEMO -- run demonstrations
%     TRIAXIAL.DOC -- this documentation
%     TRIAXIAL.TESTS -- the test suite
%     CARTPROJ -- plot a geodesic
%     HORIZON -- return point on the horizon

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  help triaxial/doc
end
