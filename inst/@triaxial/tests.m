function tests
%TRIAXIAL.TESTS  self tests for the TRIAXIAL class
%
%   TRIAXIAL.TESTS
%
%   runs a variety of tests for the TRIAIAL class.  It produces no output it
%   they are successful.

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  n = 0;
  tol = eps^(3/4);
  t = triaxial(sqrt([2,1,1/2]));
  % Generate a representative set of cartesian coordinates
  r = reshape(magic(9), 27,3);
  r = (r - mean(r(:))) / 20 .* t.axes;
  [r2, h] = t.carttocart2(r); rn = t.cart2tocart(r2, h);
  n = n + assertEquals(rn, r, tol);
  r2n = t.cart2norm(r2);
  n = n + assertEquals(r2n, r2, tol);
  ellip = t.carttoellip(r); rn = t.elliptocart(ellip);
  n = n + assertEquals(rn, r, tol);
  geod = t.carttogeod(r); rn = t.geodtocart(geod);
  n = n + assertEquals(rn, r, tol);
  param = t.carttoparam(r); rn = t.paramtocart(param);
  n = n + assertEquals(rn, r, tol);
  geocen = t.carttogeocen(r); rn = t.geocentocart(geocen);
  n = n + assertEquals(rn, r, tol);
  ellip = t.cart2toellip(r2); r2n = t.elliptocart(ellip);
  n = n + assertEquals(r2n, r2, tol);
  geod = t.cart2togeod(r2); r2n = t.geodtocart(geod);
  param = t.cart2toparam(r2); r2n = t.paramtocart(param);
  n = n + assertEquals(r2n, r2, tol);
  geocen = t.cart2togeocen(r2); r2n = t.geocentocart(geocen);
  n = n + assertEquals(r2n, r2, tol);
  rr = t.convert(r, 'cartesian', 'ellipsoidal');
  rr = t.convert(rr, 'ellipsoidal', 'geodetic');
  rr = t.convert(rr, 'geodetic', 'parametric');
  rr = t.convert(rr, 'parametric', 'geocentric');
  rr = t.convert(rr, 'geocentric', 'ellipsoidal');
  rr = t.convert(rr, 'ellipsoidal', 'cartesian');
  n = n + assertEquals(rr, r, tol);
  rr = t.convert(r2, 'cartesian2', 'ellipsoidal');
  rr = t.convert(rr, 'ellipsoidal', 'geodetic');
  rr = t.convert(rr, 'geodetic', 'parametric');
  rr = t.convert(rr, 'parametric', 'geocentric');
  rr = t.convert(rr, 'geocentric', 'ellipsoidal');
  rr = t.convert(rr, 'ellipsoidal', 'cartesian2');
  n = n + assertEquals(rr, r2, tol);

  v2 = reshape(r2, 3, 27)';
  [r2, v2] = t.cart2norm(r2, v2);
  [ellip, alp] = t.cart2toellip(r2, v2);
  [r2n, v2n] = t.elliptocart2(ellip, alp);
  n = n + assertEquals(r2n, r2, tol);
  n = n + assertEquals(v2n, v2, tol);
  assert(n == 0);
end

function n = assertEquals(u, v, d)
  n = abs(u - v) <= d;
  n = sum(~n(:));
end
