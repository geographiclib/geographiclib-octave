function [r2f, r2b] = cartproj(t, r, viewpt, varargin)
%CARTPROJ  Project and plot a curve on the surface of the ellipsoid
%
%   r2f = CARTPROJ(t, r, viewpt)
%   [r2f, r2b] = CARTPROJ(t, r, viewpt)
%   r3 = CARTPROJ(t, r, [])
%   CARTPROJ(t, r, viewpt, ...)
%   CARTPROJ(t, r, [], ...)
%
%   Input:
%     t the triaxial ellipsoid object
%     r an n x 3 array of 3d cartesian points defining the curve
%     viewpt the geodetic coordinates defining the viewing direction
%     ... additional arguments for plotting the curve
%   Output:
%     r2f an n x 2 array of projected points facing toward the viewpt
%     r2b an n x 2 array of projected points facing away from the viewpt
%     r an n x 3 array of 3d points
%
%   This routine can be used to plot geodesics produced by RECKON.  The input r
%   can contained nans; this can be used to plot disconnected lines as a single
%   array (e.g., the latitude-longitude graticule or several geodesics
%   emanating from a single point).
%
%   The viewpt is defined by the geodetic coordinates [phi, lam].  The points
%   r2f (front) and r2b (back) are in an orthographic projection from a distant
%   point above this viewpt.  The hidden points are replaced by nans.  A little
%   "fuzz" is included in the determination hidden points -- points slightly
%   over the horizon are counted as visible.  No attempt is made to interpolate
%   the point on the horizon, so the points need to be sufficiently close so
%   that the curves extend to the horizon.  There's little cost in this
%   requirement.  RECKON can return a dense curve at minimal cost and the
%   subsequent call to LINESIMP will cull unneeded points.
%
%   If viewpt is set to [], then no projection is performed and the 3d
%   points are simplified and returned.
%
%   In all cases the points (either the 2d projected points or the original 3d
%   points) are passed through LINESIMP.  This always collapses runs of
%   consecutive nans.  In addition, it applies the Ramer-Douglas-Peucker
%   algorithms to simplify the line with tolerance set to t.b * t.linesimptol.
%   If t.linesimptol < 0, then no line simplification is done (except for
%   collapsing runs of nans).
%
%   If there are no output arguments, the points are passed either to plot3
%   (if viewpt == []) or plot (otherwise), with the extra arguments passed
%   to the plotting function, e.g., to set line width or color.
%
%   See also LINESIMP, HORIZON, RECKON

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  threed = isempty(viewpt);
  plotp = nargout == 0;
  if threed
    r2f = linesimp(r, t.linesimptol, t.b * [1,1,1]);
    if plotp, plot3(r2f(:,1), r2f(:,2), r2f(:,3), varargin{:}); end
  else
    M = rotm(viewpt);
    upv = M(3,:);
    M = M(1:2,:)';
    up = r ./ t.axes.^2;
    hidden = up * upv' < -0.005;
    r2f = r;
    r2f(hidden,:) = nan;
    r2f = linesimp(r2f * M, t.linesimptol, t.b * [1,1]);
    if plotp
      plot(r2f(:,1), r2f(:,2), varargin{:})
    elseif nargout == 2
      hidden = up * upv' > 0.005;
      r2b = r;
      r2b(hidden,:) = nan;
      r2b = linesimp(r2b * M, t.linesimptol, t.b * [1,1]);
    end
  end
end
