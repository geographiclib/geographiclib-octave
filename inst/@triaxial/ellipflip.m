function [ellipn, alpn] = ellipflip(ellip, alp)
%TRIAXIAL.ELLIPFLIP  switch ellipsoidal coordinates to other sheet
%
%   ellip = TRIAXIAL.ELLIPFLIP(ellip)
%   [ellip, alp] = TRIAXIAL.ELLIPFLIP(ellip, alp)
%   [ellip, alp, switch] = TRIAXIAL.ELLIPNORM(ellip, alp, alt)
%
%   Input:
%     ellip an n x 2 array of ellipsoidal coordinates [bet, omg]
%     alp an n x 1 array of azimuths alpha
%   Output:
%     ellipn an n x 2 array of ellipsoidal coordinates [bet, omg]
%     alpn an n x 1 array of azimuths alpha.
%     flip an n x 1 logical array specifies if the sheet was switched
%
%   bet, omg, and alp are measured in degrees.  This reduces bet to the
%   range [-90, 90] and omg and alp to the range [-180, 180].  If alt is
%   set, the ranges are instead: bet in [-180, 180], omg in [0, 180].

% Copyright (c) Charles Karney (2024) <karney@alum.mit.edu>.
 
  bet = ellip(:,1);
  omg = ellip(:,2);
  betn = remx( signx(bet)*180-bet, 360, true);
  omgn = remx(               -omg, 360, true);
  ellipn = [betn, omgn];
  if nargin > 1
    alpn = remx(-signx(alp)*180+alp, 360, true);
  end

end

