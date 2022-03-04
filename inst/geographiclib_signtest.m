function geographiclib_signtest
%SIGN_TEST   The test suite for the geographiclib package
%
%   SIGN_TEST
%
%   runs a variety of tests checking the treatment of +/-0 and +/-180

  persistent octavep
  if isempty(octavep)
    octavep = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  end

  n = 0;

  i = checkAngRound(-eps/32, -eps/32);
  if i, n=n+1; fprintf('AngRound(-eps/32) fail\n'); end
  i = checkAngRound(-eps/64, -0.0   );
  if i, n=n+1; fprintf('AngRound(-eps/64) fail\n'); end
  i = checkAngRound(-  0.0 , -0.0   );
  if i, n=n+1; fprintf('AngRound(-  0.0 ) fail\n'); end
  i = checkAngRound(   0.0 , +0.0   );
  if i, n=n+1; fprintf('AngRound(   0.0 ) fail\n'); end
  i = checkAngRound( eps/64, +0.0   );
  if i, n=n+1; fprintf('AngRound( eps/64) fail\n'); end
  i = checkAngRound( eps/32, +eps/32);
  if i, n=n+1; fprintf('AngRound( eps/32) fail\n'); end
  i = checkAngRound((1-2*eps)/64, (1-2*eps)/64);
  if i, n=n+1; fprintf('AngRound((1-2*eps)/64) fail\n'); end
  i = checkAngRound((1-eps  )/64,  1.0     /64);
  if i, n=n+1; fprintf('AngRound((1-eps  )/64) fail\n'); end
  i = checkAngRound((1-eps/2)/64,  1.0     /64);
  if i, n=n+1; fprintf('AngRound((1-eps/2)/64) fail\n'); end
  i = checkAngRound((1-eps/4)/64,  1.0     /64);
  if i, n=n+1; fprintf('AngRound((1-eps/4)/64) fail\n'); end
  i = checkAngRound( 1.0     /64,  1.0     /64);
  if i, n=n+1; fprintf('AngRound( 1.0     /64) fail\n'); end
  i = checkAngRound((1+eps/2)/64,  1.0     /64);
  if i, n=n+1; fprintf('AngRound((1+eps/2)/64) fail\n'); end
  i = checkAngRound((1+eps  )/64,  1.0     /64);
  if i, n=n+1; fprintf('AngRound((1+eps  )/64) fail\n'); end
  i = checkAngRound((1+2*eps)/64, (1+2*eps)/64);
  if i, n=n+1; fprintf('AngRound((1+2*eps)/64) fail\n'); end
  i = checkAngRound((1-eps  )/32, (1-eps  )/32);
  if i, n=n+1; fprintf('AngRound((1-eps  )/32) fail\n'); end
  i = checkAngRound((1-eps/2)/32,  1.0     /32);
  if i, n=n+1; fprintf('AngRound((1-eps/2)/32) fail\n'); end
  i = checkAngRound((1-eps/4)/32,  1.0     /32);
  if i, n=n+1; fprintf('AngRound((1-eps/4)/32) fail\n'); end
  i = checkAngRound( 1.0     /32,  1.0     /32);
  if i, n=n+1; fprintf('AngRound( 1.0     /32) fail\n'); end
  i = checkAngRound((1+eps/2)/32,  1.0     /32);
  if i, n=n+1; fprintf('AngRound((1+eps/2)/32) fail\n'); end
  i = checkAngRound((1+eps  )/32, (1+eps  )/32);
  if i, n=n+1; fprintf('AngRound((1+eps  )/32) fail\n'); end
  i = checkAngRound((1-eps  )/16, (1-eps  )/16);
  if i, n=n+1; fprintf('AngRound((1-eps  )/16) fail\n'); end
  i = checkAngRound((1-eps/2)/16, (1-eps/2)/16);
  if i, n=n+1; fprintf('AngRound((1-eps/2)/16) fail\n'); end
  i = checkAngRound((1-eps/4)/16,  1.0     /16);
  if i, n=n+1; fprintf('AngRound((1-eps/4)/16) fail\n'); end
  i = checkAngRound( 1.0     /16,  1.0     /16);
  if i, n=n+1; fprintf('AngRound( 1.0     /16) fail\n'); end
  i = checkAngRound((1+eps/4)/16,  1.0     /16);
  if i, n=n+1; fprintf('AngRound((1+eps/4)/16) fail\n'); end
  i = checkAngRound((1+eps/2)/16,  1.0     /16);
  if i, n=n+1; fprintf('AngRound((1+eps/2)/16) fail\n'); end
  i = checkAngRound((1+eps  )/16, (1+eps  )/16);
  if i, n=n+1; fprintf('AngRound((1+eps  )/16) fail\n'); end
  i = checkAngRound((1-eps  )/ 8, (1-eps  )/ 8);
  if i, n=n+1; fprintf('AngRound((1-eps  )/ 8) fail\n'); end
  i = checkAngRound((1-eps/2)/ 8, (1-eps/2)/ 8);
  if i, n=n+1; fprintf('AngRound((1-eps/2)/ 8) fail\n'); end
  i = checkAngRound((1-eps/4)/ 8,  1.0     / 8);
  if i, n=n+1; fprintf('AngRound((1-eps/4)/ 8) fail\n'); end
  i = checkAngRound((1+eps/2)/ 8,  1.0     / 8);
  if i, n=n+1; fprintf('AngRound((1+eps/2)/ 8) fail\n'); end
  i = checkAngRound((1+eps  )/ 8, (1+eps  )/ 8);
  if i, n=n+1; fprintf('AngRound((1+eps  )/ 8) fail\n'); end
  i = checkAngRound( 1-eps      ,  1-eps      );
  if i, n=n+1; fprintf('AngRound( 1-eps      ) fail\n'); end
  i = checkAngRound( 1-eps/2    ,  1-eps/2    );
  if i, n=n+1; fprintf('AngRound( 1-eps/2    ) fail\n'); end
  i = checkAngRound( 1-eps/4    ,  1          );
  if i, n=n+1; fprintf('AngRound( 1-eps/4    ) fail\n'); end
  i = checkAngRound( 1.0        ,  1          );
  if i, n=n+1; fprintf('AngRound( 1.0        ) fail\n'); end
  i = checkAngRound( 1+eps/4    ,  1          );
  if i, n=n+1; fprintf('AngRound( 1+eps/4    ) fail\n'); end
  i = checkAngRound( 1+eps/2    ,  1          );
  if i, n=n+1; fprintf('AngRound( 1+eps/2    ) fail\n'); end
  i = checkAngRound( 1+eps      ,  1+  eps    );
  if i, n=n+1; fprintf('AngRound( 1+eps      ) fail\n'); end
  i = checkAngRound( 90.0-64*eps,  90-64*eps  );
  if i, n=n+1; fprintf('AngRound( 90.0-64*eps) fail\n'); end
  i = checkAngRound( 90.0-32*eps,  90         );
  if i, n=n+1; fprintf('AngRound( 90.0-32*eps) fail\n'); end
  i = checkAngRound( 90.0       ,  90         );
  if i, n=n+1; fprintf('AngRound( 90.0       ) fail\n'); end

  i = checksincosd(-  inf,  nan,  nan);
  if i, n=n+1; fprintf('sincosd(-  inf) fail\n'); end
  i = checksincosd(-810.0, -1.0, +0.0);
  if i, n=n+1; fprintf('sincosd(-810.0) fail\n'); end
  i = checksincosd(-720.0, -0.0, +1.0);
  if i, n=n+1; fprintf('sincosd(-720.0) fail\n'); end
  i = checksincosd(-630.0, +1.0, +0.0);
  if i, n=n+1; fprintf('sincosd(-630.0) fail\n'); end
  i = checksincosd(-540.0, -0.0, -1.0);
  if i, n=n+1; fprintf('sincosd(-540.0) fail\n'); end
  i = checksincosd(-450.0, -1.0, +0.0);
  if i, n=n+1; fprintf('sincosd(-450.0) fail\n'); end
  i = checksincosd(-360.0, -0.0, +1.0);
  if i, n=n+1; fprintf('sincosd(-360.0) fail\n'); end
  i = checksincosd(-270.0, +1.0, +0.0);
  if i, n=n+1; fprintf('sincosd(-270.0) fail\n'); end
  i = checksincosd(-180.0, -0.0, -1.0);
  if i, n=n+1; fprintf('sincosd(-180.0) fail\n'); end
  i = checksincosd(- 90.0, -1.0, +0.0);
  if i, n=n+1; fprintf('sincosd(- 90.0) fail\n'); end
  i = checksincosd(-  0.0, -0.0, +1.0);
  if i, n=n+1; fprintf('sincosd(-  0.0) fail\n'); end
  i = checksincosd(+  0.0, +0.0, +1.0);
  if i, n=n+1; fprintf('sincosd(+  0.0) fail\n'); end
  i = checksincosd(+ 90.0, +1.0, +0.0);
  if i, n=n+1; fprintf('sincosd(+ 90.0) fail\n'); end
  i = checksincosd(+180.0, +0.0, -1.0);
  if i, n=n+1; fprintf('sincosd(+180.0) fail\n'); end
  i = checksincosd(+270.0, -1.0, +0.0);
  if i, n=n+1; fprintf('sincosd(+270.0) fail\n'); end
  i = checksincosd(+360.0, +0.0, +1.0);
  if i, n=n+1; fprintf('sincosd(+360.0) fail\n'); end
  i = checksincosd(+450.0, +1.0, +0.0);
  if i, n=n+1; fprintf('sincosd(+450.0) fail\n'); end
  i = checksincosd(+540.0, +0.0, -1.0);
  if i, n=n+1; fprintf('sincosd(+540.0) fail\n'); end
  i = checksincosd(+630.0, -1.0, +0.0);
  if i, n=n+1; fprintf('sincosd(+630.0) fail\n'); end
  i = checksincosd(+720.0, +0.0, +1.0);
  if i, n=n+1; fprintf('sincosd(+720.0) fail\n'); end
  i = checksincosd(+810.0, +1.0, +0.0);
  if i, n=n+1; fprintf('sincosd(+810.0) fail\n'); end
  i = checksincosd(+  inf,  nan,  nan);
  if i, n=n+1; fprintf('sincosd(+  inf) fail\n'); end
  i = checksincosd(   nan,  nan,  nan);
  if i, n=n+1; fprintf('sincosd(   nan) fail\n'); end

  [s1, c1] = sincosdx(9);
  [s2, c2] = sincosdx(81);
  [s3, c3] = sincosdx(-123456789);
  i = equiv(s1, c2) + equiv(s1, s3) + equiv(c1, s2) + equiv(c1, -c3);
  if i, n=n+1; fprintf('sincosd accuracy fail\n'); end

  i = checkatan2d(+0.0 , +0.0 , +0.0 );
  if i, n=n+1; fprintf('atan2d(+0.0 , +0.0 ) fail\n'); end
  if octavep
    % MATLAB gets these wrong, but we don't care about the 0,0 pairs
    i = checkatan2d(+0.0 , -0.0 , +180 );
    if i, n=n+1; fprintf('atan2d(+0.0 , -0.0 ) fail\n'); end
    i = checkatan2d(-0.0 , -0.0 , -180 );
    if i, n=n+1; fprintf('atan2d(-0.0 , -0.0 ) fail\n'); end
    i = checkatan2d(-0.0 , +0.0 , -0.0 );
    if i, n=n+1; fprintf('atan2d(-0.0 , +0.0 ) fail\n'); end
  end
  i = checkatan2d(+0.0 , -1.0 , +180 );
  if i, n=n+1; fprintf('atan2d(+0.0 , -1.0 ) fail\n'); end
  i = checkatan2d(-0.0 , -1.0 , -180 );
  if i, n=n+1; fprintf('atan2d(-0.0 , -1.0 ) fail\n'); end
  i = checkatan2d(+0.0 , +1.0 , +0.0 );
  if i, n=n+1; fprintf('atan2d(+0.0 , +1.0 ) fail\n'); end
  i = checkatan2d(-0.0 , +1.0 , -0.0 );
  if i, n=n+1; fprintf('atan2d(-0.0 , +1.0 ) fail\n'); end
  i = checkatan2d(-1.0 , +0.0 ,  -90 );
  if i, n=n+1; fprintf('atan2d(-1.0 , +0.0 ) fail\n'); end
  i = checkatan2d(-1.0 , -0.0 ,  -90 );
  if i, n=n+1; fprintf('atan2d(-1.0 , -0.0 ) fail\n'); end
  i = checkatan2d(+1.0 , +0.0 ,  +90 );
  if i, n=n+1; fprintf('atan2d(+1.0 , +0.0 ) fail\n'); end
  i = checkatan2d(+1.0 , -0.0 ,  +90 );
  if i, n=n+1; fprintf('atan2d(+1.0 , -0.0 ) fail\n'); end
  i = checkatan2d(+1.0 ,  -inf, +180 );
  if i, n=n+1; fprintf('atan2d(+1.0 ,  -inf) fail\n'); end
  i = checkatan2d(-1.0 ,  -inf, -180 );
  if i, n=n+1; fprintf('atan2d(-1.0 ,  -inf) fail\n'); end
  i = checkatan2d(+1.0 ,  +inf, +0.0 );
  if i, n=n+1; fprintf('atan2d(+1.0 ,  +inf) fail\n'); end
  i = checkatan2d(-1.0 ,  +inf, -0.0 );
  if i, n=n+1; fprintf('atan2d(-1.0 ,  +inf) fail\n'); end
  i = checkatan2d( +inf, +1.0 ,  +90 );
  if i, n=n+1; fprintf('atan2d( +inf, +1.0 ) fail\n'); end
  i = checkatan2d( +inf, -1.0 ,  +90 );
  if i, n=n+1; fprintf('atan2d( +inf, -1.0 ) fail\n'); end
  i = checkatan2d( -inf, +1.0 ,  -90 );
  if i, n=n+1; fprintf('atan2d( -inf, +1.0 ) fail\n'); end
  i = checkatan2d( -inf, -1.0 ,  -90 );
  if i, n=n+1; fprintf('atan2d( -inf, -1.0 ) fail\n'); end
  i = checkatan2d( +inf,  -inf, +135 );
  if i, n=n+1; fprintf('atan2d( +inf,  -inf) fail\n'); end
  i = checkatan2d( -inf,  -inf, -135 );
  if i, n=n+1; fprintf('atan2d( -inf,  -inf) fail\n'); end
  i = checkatan2d( +inf,  +inf,  +45 );
  if i, n=n+1; fprintf('atan2d( +inf,  +inf) fail\n'); end
  i = checkatan2d( -inf,  +inf,  -45 );
  if i, n=n+1; fprintf('atan2d( -inf,  +inf) fail\n'); end
  i = checkatan2d(  nan, +1.0 ,  nan );
  if i, n=n+1; fprintf('atan2d(  nan, +1.0 ) fail\n'); end
  i = checkatan2d(+1.0 ,   nan,  nan );
  if i, n=n+1; fprintf('atan2d(+1.0 ,   nan) fail\n'); end

  s = 7e-16; i = equiv( atan2dx(s, -1), 180 - atan2dx(s, 1) );
  if i, n=n+1; fprintf('atan2d accuracy fail\n'); end

  i = checksum(+9.0, -9.0, +0.0 );
  if i, n=n+1; fprintf('sum(+9.0, -9.0) fail\n'); end
  i = checksum(-9.0, +9.0, +0.0 );
  if i, n=n+1; fprintf('sum(-9.0, +9.0) fail\n'); end
  i = checksum(-0.0, +0.0, +0.0 );
  if i, n=n+1; fprintf('sum(-0.0, +0.0) fail\n'); end
  i = checksum(+0.0, -0.0, +0.0 );
  if i, n=n+1; fprintf('sum(+0.0, -0.0) fail\n'); end
  i = checksum(-0.0, -0.0, -0.0 );
  if i, n=n+1; fprintf('sum(-0.0, -0.0) fail\n'); end
  i = checksum(+0.0, +0.0, +0.0 );
  if i, n=n+1; fprintf('sum(+0.0, +0.0) fail\n'); end

  i = checkAngNormalize(-900.0, -180 );
  if i, n=n+1; fprintf('AngNormalize(-900.0) fail\n'); end
  i = checkAngNormalize(-720.0, -0.0 );
  if i, n=n+1; fprintf('AngNormalize(-720.0) fail\n'); end
  i = checkAngNormalize(-540.0, -180 );
  if i, n=n+1; fprintf('AngNormalize(-540.0) fail\n'); end
  i = checkAngNormalize(-360.0, -0.0 );
  if i, n=n+1; fprintf('AngNormalize(-360.0) fail\n'); end
  i = checkAngNormalize(-180.0, -180 );
  if i, n=n+1; fprintf('AngNormalize(-180.0) fail\n'); end
  i = checkAngNormalize(  -0.0, -0.0 );
  if i, n=n+1; fprintf('AngNormalize(  -0.0) fail\n'); end
  i = checkAngNormalize(  +0.0, +0.0 );
  if i, n=n+1; fprintf('AngNormalize(  +0.0) fail\n'); end
  i = checkAngNormalize( 180.0, +180 );
  if i, n=n+1; fprintf('AngNormalize( 180.0) fail\n'); end
  i = checkAngNormalize( 360.0, +0.0 );
  if i, n=n+1; fprintf('AngNormalize( 360.0) fail\n'); end
  i = checkAngNormalize( 540.0, +180 );
  if i, n=n+1; fprintf('AngNormalize( 540.0) fail\n'); end
  i = checkAngNormalize( 720.0, +0.0 );
  if i, n=n+1; fprintf('AngNormalize( 720.0) fail\n'); end
  i = checkAngNormalize( 900.0, +180 );
  if i, n=n+1; fprintf('AngNormalize( 900.0) fail\n'); end

  i = checkAngDiff(+  0.0, +  0.0,   +0.0 );
  if i, n=n+1; fprintf('AngDiff(+  0.0, +  0.0) fail\n'); end
  i = checkAngDiff(+  0.0, -  0.0,   -0.0 );
  if i, n=n+1; fprintf('AngDiff(+  0.0, -  0.0) fail\n'); end
  i = checkAngDiff(-  0.0, +  0.0,   +0.0 );
  if i, n=n+1; fprintf('AngDiff(-  0.0, +  0.0) fail\n'); end
  i = checkAngDiff(-  0.0, -  0.0,   +0.0 );
  if i, n=n+1; fprintf('AngDiff(-  0.0, -  0.0) fail\n'); end
  i = checkAngDiff(+  5.0, +365.0,   +0.0 );
  if i, n=n+1; fprintf('AngDiff(+  5.0, +365.0) fail\n'); end
  i = checkAngDiff(+365.0, +  5.0,   -0.0 );
  if i, n=n+1; fprintf('AngDiff(+365.0, +  5.0) fail\n'); end
  i = checkAngDiff(+  5.0, +185.0, +180.0 );
  if i, n=n+1; fprintf('AngDiff(+  5.0, +185.0) fail\n'); end
  i = checkAngDiff(+185.0, +  5.0, -180.0 );
  if i, n=n+1; fprintf('AngDiff(+185.0, +  5.0) fail\n'); end
  i = checkAngDiff( +eps , +180.0, +180.0 );
  if i, n=n+1; fprintf('AngDiff( +eps , +180.0) fail\n'); end
  i = checkAngDiff( -eps , +180.0, -180.0 );
  if i, n=n+1; fprintf('AngDiff( -eps , +180.0) fail\n'); end
  i = checkAngDiff( +eps , -180.0, +180.0 );
  if i, n=n+1; fprintf('AngDiff( +eps , -180.0) fail\n'); end
  i = checkAngDiff( -eps , -180.0, -180.0 );
  if i, n=n+1; fprintf('AngDiff( -eps , -180.0) fail\n'); end

  % azimuth of geodesic line with points on equator determined by signs of
  % latitude
  % lat1 lat2 azi1/2
  C = [ +0, -0, 180;
        -0, +0, 0 ];
  [~, azi1, azi2] = geoddistance(C(:,1), 0, C(:,2), 0);
  i = equiv(azi1, C(:,3)) + equiv(azi2, C(:,3));
  if i, n=n+1; fprintf('coincident points on equator fail\n'); end

  % Does the nearly antipodal equatorial solution go north or south?
  % lat1 lat2 azi1 azi2
  C = [ +0, +0,  56, 124;
        -0, -0, 124,  56];
  [~, azi1, azi2] = geoddistance(C(:,1), 0, C(:,2), 179.5);
  i = assertEquals(azi1, C(:,3), 1) + assertEquals(azi2, C(:,4), 1);
  if i, n=n+1; fprintf('nearly antipodal points on equator fail\n'); end

  % How does the exact antipodal equatorial path go N/S + E/W
  % lat1 lat2 lon2 azi1 azi2
  C = [ +0, +0, +180,   +0, +180;
        -0, -0, +180, +180,   +0;
        +0, +0, -180,   -0, -180;
        -0, -0, -180, -180,   -0 ];
  [~, azi1, azi2] = geoddistance(C(:,1), 0, C(:,2), C(:,3));
  i = equiv(azi1, C(:,4)) + equiv(azi2, C(:,5));
  if i, n=n+1; fprintf('antipodal points on equator fail\n'); end

  % Antipodal points on the equator with prolate ellipsoid
  % lon2 azi1/2
  C = [+180, +90;
       -180, -90];
  [~, azi1, azi2] = geoddistance(0, 0, 0, C(:,1), [6.4e6, flat2ecc(-1/300)]);
  i = equiv(azi1, C(:,2)) + equiv(azi2, C(:,2));
  if i, n=n+1; fprintf('`antipodal points on equator, prolate, fail\n'); end

  % azimuths = +/-0 and +/-180 for the direct problem
  % azi1, lon2, azi2
  C = [  +0, +180, +180;
         -0, -180, -180;
       +180, +180,   +0;
       -180, -180,   -0 ];
  [~, lon2, azi2] = geodreckon(0, 0, 15e6, C(:,1), 2);
  i = equiv(lon2, C(:,2)) + equiv(azi2, C(:,3));
  if i, n=n+1; fprintf('geodreckon azi1 = +/-0 +/-180, fail\n'); end

  % lat = +/-0 in utmups_fwd
  % lat y northp
  C = [ +0,  0,   1;
        -0, 10e6, 0 ];
  [x, y, zone, northp] = utmups_fwd(C(:,1), 3);
  i = equiv(y, C(:,2)) + equiv(northp, C(:,3));
  if i, n=n+1; fprintf('utmups_fwd lat = +/-0, fail\n'); end
  mgrs = mgrs_fwd(x, y, zone, northp, 2);
  i = ~strcmp(mgrs{1,1}, '31NEA0000') + ~strcmp(mgrs{2,1}, '31MEV0099');
  if i, n=n+1; fprintf('mgrs_fwd lat = +/-0, fail\n'); end

  assert(n == 0);
end

function n = assertEquals(x, y, d)
  n = abs(x - y) <= d;
  n = sum(~n(:));
end

function n = equiv(x, y)
  n = (isnan(x) & isnan(y)) | (x == y & signbitx(x) == signbitx(y));
  n = sum(~n(:));
end

function n = checkAngRound(x, y)
  z = AngRound(x);
  n = equiv(z, y);
end

function n = checksincosd(x, s, c)
  [ss, cc] = sincosdx(x);
  n = equiv(ss, s) + equiv(cc, c);
end

function n = checkatan2d(y, x, a)
  aa = atan2dx(y, x);
  n = equiv(aa, a);
end

function n = checksum(x, y, s)
  ss = sumx(x, y);
  n = equiv(ss, s);
end

function n = checkAngNormalize(x, y)
  z = AngNormalize(x);
  n = equiv(z, y);
end

function n = checkAngDiff(x, y, d)
  dd = AngDiff(x, y);
  n = equiv(dd, d);
end

%!test geographiclib_signtest
