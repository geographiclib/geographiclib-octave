function ptsout = linesimp(ptsin, tol, scal)
%LINESIMP  simplify a line using the Ramer-Douglas-Peucker method
%
%   ptsout = LINESIMP(ptsin, tol)
%   ptsout = LINESIMP(ptsin, tol, scal)
%
%   Input:
%     ptsin an n x 2 or n x 3 array of input points
%     tol the tolerane used to determine whether points can be removed
%     scal a vector of length 2 o3 3 used to prescale the points
%   Output:
%     ptsout an m x 2 or m x 3 array of input points
%
%   Simplify a polyline by removing redundant points.  This routine can deal
%   with points in 2d or 3d.  The routine is nan-aware; consecutive runs of
%   nans are collapsed to a single entry and initial and final runs of nans
%   are removed.  Particularly for 2d data, it's useful to be able to scale
%   the data before the algorithm is applied using scale.  For example:
%
%     x = [0:360];;
%     xyin = [x,sind(x)];
%     xyout = linesimp(xyin, 0.001, [180/pi, 1]);
%
%   If scal is omitted, it is assumed to be [1,1] or [1,1,1].  If tol < 0, it
%   just returns the original points with nans collapsed.

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

% Matlab central contains several package which implement the same
% algorithm.  They all suffer from at least one of the following defects:
%
%   * nans are not handled
%   * the criterion for removing points may remove significant points beyond
%     the endpoints
%   * no prescaling of the date is offered
%
% As of 2024-03-29...
%
% Doesn't handle nans:
%   https://www.mathworks.com/matlabcentral/fileexchange/101455
%
% Doesn't correctly compute distance beyond endpoints:
%   https://www.mathworks.com/matlabcentral/fileexchange/21132
%
% Doesn't handle nans + doesn't correctly compute distance beyond endpoints:
%   https://www.mathworks.com/matlabcentral/fileexchange/36229
%   https://www.mathworks.com/matlabcentral/fileexchange/41986
%   https://www.mathworks.com/matlabcentral/fileexchange/61046
%
% No code:
%   https://www.mathworks.com/matlabcentral/fileexchange/72452

% first break into non-nan segments
  dim = size(ptsin, 2);
  if nargin < 3, scal = ones(1, dim); end
  scal = scal(:)';
  ptsin = ptsin ./ scal;
  b = diff([1;isnan(sum(ptsin,2));1]);
  inds = find(abs(b) == 1);
  nanpt = nan(1, dim);
  ptsout = zeros(0, dim);
  for j = 1:length(inds)/2
    if j > 1
      ptsout = [ptsout; nanpt];         %#ok<*AGROW>
    end
    ptsout = [ptsout; simp1(ptsin(inds(2*j-1):inds(2*j)-1,:), tol)];
  end
  ptsout = ptsout .* scal;
end

function ptsout = simp1(ptsin, tol)
% no nans in pts
  if tol < 0
    ptsout = ptsin;
    return
  end
  n = size(ptsin, 1);
  if n <= 2
    if n == 2 && all(ptsin(1,:) == ptsin(n,:))
      ptsout = ptsin(1,:);
    else
      ptsout = ptsin;
    end
    return
  end
  [dist, k] = maxdist(ptsin);
  if dist <= tol
    if all(ptsin(1,:) == ptsin(n,:))
      ptsout = ptsin(1,:);
    else
      ptsout = ptsin([1,n],:);
    end
    return
  end
  ptsa = simp1(ptsin(1:k,:), tol);
  ptsb = simp1(ptsin(k:n,:), tol);
  ptsout = [ptsa; ptsb(2:end,:)];
end

function [dist, k] = maxdist(pts)
  n = size(pts, 1);
  dim = size(pts, 2);
  pts = pts - pts(1,:);
  v =  pts(n,:);
  r = vecabs(v);
  if r > 0
    vn = vecunit(v);
    if dim == 2
      d = vecabs(veccross([pts, zeros(n, 1)], [vn, 0]));
    else
      d = vecabs(veccross(pts, vn));
    end
    x = vecdot(pts, vn);
    l = x > r/2;
    pts(l,:) = r * vn - pts(l,:);
    l = vecdot(vn, pts) < 0;
    d(l) = vecabs(pts(l,:));
  else
    d = vecabs(pts);
  end
  [dist, k] = max(d);
end

function y = vecabs(x)
% Input:
%  x = n x 3 vectors
% Output:
%  y = n x 1 of absolute values
%
% Equivalent to vecnorm(x,2,2)
  y = sqrt(sum(x.^2, 2));
end

function z = vecdot(x, y)
% Input:
%  x,y = n x 3 vectors
% Output:
%  z = n x 1 of dot products
%
% Like dot(x,y,2) but unlike the builtin routine, this can deal with N . 1
% or 1 . N cases.
  z = sum(x .* y, 2);
end

function y = vecunit(x)
% Input:
%  x = n x 3 vectors
% Output:
%  y = n x 3 vectors converted to unit vectors
  z = vecabs(x);
  z(z == 0) = 1;
  y = x ./ z;
end

function z = veccross(x,y)
% Input:
%  x,y = n x 3 vectors
% Output:
%  z = n x 3 vector of cross products
%
% Like cross(x,y,2), but unlike the builtin routine, this can deal with n x
% 1 or 1 x n cases.
  z = [x(:,2).*y(:,3) - x(:,3).*y(:,2), ...
       x(:,3).*y(:,1) - x(:,1).*y(:,3), ...
       x(:,1).*y(:,2) - x(:,2).*y(:,1)];
end
