function [ellipn, alpn, flip] = ellipnorm(ellip, alp, alt)
%TRIAXIAL.ELLIPNORM  reduce ellipsoidal coordinates to convention range
%
%   ellip = TRIAXIAL.ELLIPNORM(ellip)
%   [ellip, alp] = TRIAXIAL.ELLIPNORM(ellip, alp)
%   [ellip, alp, switch] = TRIAXIAL.ELLIPNORM(ellip, alp, alt)
%
%   Input:
%     ellip an n x 2 array of ellipsoidal coordinates [bet, omg]
%     alp an n x 1 array of azimuths alpha
%     alt a flag, if true use alternate ranges for bet and omg (default false)
%   Output:
%     ellipn an n x 2 array of ellipsoidal coordinates [bet, omg]
%     alpn an n x 1 array of azimuths alpha.
%     flip an n x 1 logical array specifies if the sheet was switched
%
%   bet, omg, and alp are measured in degrees.  This reduces bet to the
%   range [-90, 90] and omg and alp to the range [-180, 180].  If alt is
%   set, the ranges are instead: bet in [-180, 180], omg in [0, 180].

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.

  if nargin < 2, alp = nan; end
  if nargin < 3, alt = 0; end
  bet = ellip(:,1);
  omg = ellip(:,2);
  bet = remx(bet, 360, true);
  omg = remx(omg, 360, true);
  alp = remx(alp, 360, true);
  if ~alt
    % Flip if abs(bet) > 90 or abs(bet) == 90 and omg < 0
    l =  abs(bet) > 90 | (abs(bet) == 90 & omg < 0);
    if nargin >= 2
      % at northern umblic, return alp in (-180,-90] + (90, 180]
      % at southern umblic, return alp in (-90, 90]
      l = l | ...
          ((omg ==   0 | omg == 180) & ...
           (bet ==  90 &  (alp > -90 & alp <= 90)) | ...
           (bet == -90 & ~(alp > -90 & alp <= 90)));
    end
  else
    % Flip if omg < 0
    l =  omg < 0;
    if nargin >= 2
      % at northern umblic, return alp in (-180,-90] + (90, 180]
      % at southern umblic, return alp in (-90, 90]
      l = l | ...
          ((omg ==   0 | omg == 180) & ...
           (bet ==  90 &  (alp > -90 & alp <= 90)) | ...
           (bet == -90 & ~(alp > -90 & alp <= 90)));
    end
  end
  ellipn = [bet, omg];
  alpn = alp;
  if nargin >= 2
    [ellipn(l, :), alpn(l)] = triaxial.ellipflip(ellipn(l, :), alpn(l));
  else
    ellipn(l, :) = triaxial.ellipflip(ellipn(l, :));
  end
  flip = l;
end
