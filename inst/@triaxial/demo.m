function demo(n)
%TRIAXIAL.DEMO  plot some geodesics using the TRIAXIAL class
%
%   TRIAXIAL.DEMO
%   TRIAXIAL.DEMO(n)
%
%   with 0 <= n < 20 runs one of several of demonsrations of the TRIAXIAL
%   class.  If n is not specified, a list of possible demonstrations is
%   displayed.  Each demo creates a plot using the current figure.  Here is
%   the approximate time to run each demo.
%
%              time (s)
%      demo MATLAB Octave
%        0    0.1    0.4
%        1    2.7  129.3
%        2    4.6  188.9
%        3    2.2   78.1
%        4    4.8  188.1
%        5    1.3   46.5
%        6    1.8   66.7
%        7    2.0   78.1
%        8    3.1  129.3
%        9    0.1    0.6
%
%   N.B. Octave is about 40 times slower than MATLAB.  This is due to the
%   different implementations of the ODE solver.

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  titles = { 'ellipsoidal graticule';
             'circumpolar geodesic a';
             'circumpolar geodesic b';
             'transpolar geodesic a';
             'transpolar geodesic b';
             'umbilical geodesic';
             'cut locus';
             'astroid';
             'geodesic circles';
             'geodetic and parametric graticules';
             };
  if nargin == 0 || n < 0 || n >= 20
    demodoc(titles);
    return;
  end
  red = [179,27,27]/255;
  blue = [0,19,66]/100;
  green = [9,55,27]/100;
  threed = floor(n/10);
  m = n - 10*threed;
  ds = 0.025;
  dang = 2;
  ang = [-180:dang:180,nan]';
  if threed
    viewpt = [];
  else
    viewpt = [40, 30];
  end
  t = triaxial([1.01, 1, 0.8]);
  switch m
    case 0
      plotstart(t, viewpt);
      betlist = [-5:-1,1:5] * 15;
      [bet, omg] = meshgrid(betlist, ang);
      t.cartproj(t.elliptocart2([bet(:), omg(:)]), viewpt, 'Color', blue);
      omglist = [1:5,7:11] * 15;
      [omg, bet] = meshgrid(omglist, ang);
      t.cartproj(t.elliptocart2([bet(:), omg(:)]), viewpt, 'Color', green);
      princell = [0*ang, ang; ang, 0*ang; ang, 0*ang+90];
      t.cartproj(t.paramtocart2(princell), viewpt, 'Color', red);
      plotend(t, viewpt);
      title(titles{m+1});
    case {1, 2, 3, 4, 5}
      dat = [ 45.10472891 0 90 321.89789588;
              87.47973017 0 90 491.63487824;
              90 39.89350071 180 232.65487487;
              90 9.96933777 180 508.78175722;
              90 0 135 142.63587003 ];
      [pos1, dir1, s12] = deal(dat(m,1:2), dat(m,3), dat(m,4));
      npts = ceil(s12/ds);
      s12 = (-npts:npts)' * s12/npts;
      pos2 = t.reckon(pos1, dir1, s12);
      plotstart(t, viewpt);
      t.cartproj(t.elliptocart2(pos2), viewpt, 'Color', blue);
      princell = [0*ang, ang; ang, 0*ang; ang, 0*ang+90];
      t.cartproj(t.paramtocart2(princell), viewpt, 'Color', red);
      plotend(t, viewpt);
      title(titles{m+1});
    case 6
      % pos1 = t.cart2toellip(-t.geodtocart(viewpt))
      ellip1 = [-34.4575, -149.9417]; bet2 = - ellip1(1);
      plotstart(t, viewpt);
      r = zeros(0, 3);
      omg2 = [];
      for alp1 = 90:5:270
        if mod(alp1, 90) == 0
          s1 = 0;
        elseif mod(alp1, 30) == 0
          s1 = t.c/8;
        elseif mod(alp1, 10) == 0
          s1 = t.c/4;
        else
          s1 = t.c/3;
        end
        [r1, v1] = t.elliptocart2(ellip1, alp1);
        if mod(alp1-90, 180) == 0
          cond = [7, 0, -1];
        else
          cond = [11, bet2, 1];
        end
        [~,~, s12] = t.hybrid(r1, v1, cond);
        npts = ceil((s12-s1)/ds);
        s12 = linspace(s1, s12, npts)';
        r2 = t.reckon(r1, v1, s12);
        r = [r; r2];                    %#ok<*AGROW>
        if mod(alp1-90, 180) ~= 0
          alp1a = 180 - alp1;
          [r1a, v1a] = t.elliptocart2(ellip1, alp1a);
          r2a = t.reckon(r1a, v1a, s12);
          r = [r; r2a];
        else
          ellip2 = t.cart2toellip(r2(end,:));
          omg2 = [omg2; ellip2(2)];
        end
      end
      omg2 = linspace(omg2(1), omg2(2), ceil(diff(omg2)/dang))';
      cut = [bet2+0*omg2, omg2];
      t.cartproj(t.elliptocart2(cut), viewpt, 'Color', green, 'LineWidth', 2);
      princell = [0*ang, ang; ang, 0*ang; ang, 0*ang+90];
      t.cartproj(r, viewpt, 'Color', blue);
      t.cartproj(t.paramtocart2(princell), viewpt, 'Color', red);
      plotend(t, viewpt);
      title(titles{m+1});
    case 7
      % like 6 but extend to opposite line of long
      % pos1 = t.cart2toellip(-t.geodtocart(viewpt))
      ellip1 = [-34.4575, -149.9417]; omg2 = 180 + ellip1(2);
      plotstart(t, viewpt);
      r = zeros(0, 3);
      bet2 = [];
      for alp1 = 0:5:180
        if mod(alp1, 90) == 0
          s1 = 0;
        elseif mod(alp1, 30) == 0
          s1 = t.c/8;
        elseif mod(alp1, 10) == 0
          s1 = t.c/4;
        else
          s1 = t.c/3;
        end
        [r1, v1] = t.elliptocart2(ellip1, alp1);
        if mod(alp1, 180) == 0
          cond = [7, 0, -1];
        else
          cond = [12, omg2, 1];
        end
        [~,~, s12] = t.hybrid(r1, v1, cond, 1);
        npts = ceil((s12-s1)/ds);
        s12 = linspace(s1, s12, npts)';
        r2 = t.reckon(r1, v1, s12);
        r = [r; r2];
        if mod(alp1, 180) ~= 0
          alp1a = -alp1;
          [r1a, v1a] = t.elliptocart2(ellip1, alp1a);
          r2a = t.reckon(r1a, v1a, s12);
          r = [r; r2a];
        else
          ellip2 = t.cart2toellip(r2(end,:));
          bet2 = [bet2; ellip2(1)];
        end
      end
      bet2 = linspace(bet2(1), bet2(2), ceil(abs(diff(bet2))/dang))';
      cut = [bet2, omg2+0*bet2];
      t.cartproj(t.elliptocart2(cut), viewpt, 'Color', green, 'LineWidth', 2);
      princell = [0*ang, ang; ang, 0*ang; ang, 0*ang+90];
      t.cartproj(r, viewpt, 'Color', blue);
      t.cartproj(t.paramtocart2(princell), viewpt, 'Color', red);
      plotend(t, viewpt);
      title(titles{m+1});
    case 8
      % pos1 = t.cart2toellip(-t.geodtocart(viewpt))
      ellip1 = [-34.4575, -149.9417];
      ellip2 = [-ellip1(1), 180 + ellip1(2)];
      plotstart(t, viewpt);
      s12 = t.distance(ellip1, ellip2);
      s12 = (5:40)/30 * s12;
      alp1 = [-180:dang:179.9,nan,nan]';
      [r1, v1] = t.elliptocart2(ellip1, alp1);
      n = size(v1, 1);
      k = length(s12);
      r2 = t.reckon(r1, v1, s12);
      % This would plot the geodesic lines.
      %   t.cartproj(r2(ind(:),:), viewpt, 'Color', red);
      % Instead, we want to plot the geodesic circles.  First close the
      % circles in the positions of the first nan in alp1.
      r2((n-2)*k+(1:k), :) = r2(1:k, :);
      % Generate permutation to access r2 to yield the circles.
      ind = reshape(1:k*n, k, n)';
      t.cartproj(r2(ind(:),:), viewpt, 'Color', red);
      plotend(t, viewpt);
      title(titles{m+1});
    case 9
      tc = triaxial(sqrt([2, 1, 1/2]));
      plotstart(tc, viewpt);
      philist = [-5:-1,1:5]* 15;
      [phi, lam] = meshgrid(philist, ang);
      tc.cartproj(tc.geodtocart([phi(:), lam(:)]), viewpt, 'Color', blue);
      tc.cartproj(tc.paramtocart2([phi(:), lam(:)]), viewpt, 'Color', green);
      lamlist = [-11:-7,-5:-1,1:5,7:11]*15;
      anga = [-75:2:75,nan]';
      [lam, phi] = meshgrid(lamlist, anga);
      tc.cartproj(tc.geodtocart([phi(:), lam(:)]), viewpt, 'Color', blue);
      tc.cartproj(tc.paramtocart2([phi(:), lam(:)]), viewpt, 'Color', blue);
      princell = [0*ang, ang; ang, 0*ang; ang, 0*ang+90];
      tc.cartproj(tc.paramtocart2(princell), viewpt, 'Color', red);
      plotend(tc, viewpt);
      title(titles{m+1});
    otherwise
      demodoc(titles);
      error('Unknown demo');
  end
end

function plotstart(t, viewpt)
  if isempty(viewpt)
    a = t.axes - 0.001*t.c;
    [X,Y,Z] = ellipsoid(0,0,0, a(1), a(2), a(3), 30);
    C = zeros(size(X,1), size(X,2), 3) + 0.85;
    surf(X, Y, Z, C, 'FaceAlpha', 0.85);
  else
    ang = [-180:2:180,nan]';
    t.cartproj(t.horizon(ang, viewpt), viewpt, 'k');
  end
  hold on;
end

function plotend(~, viewpt)
  if isempty(viewpt)
    shading flat; axis vis3d;
  end
  axis off; axis equal; hold off;
end

function demodoc(titles)
  disp('call triaxial.demo(n) with n in [0:19];');
  disp('for n = 0:9, 2d plots:');
  for i = 0:9
    disp(['  ', int2str(i), ' display ', titles{i+1}]);
  end
  disp('for n = 10:19, the corresponding 3d plots.');
end
