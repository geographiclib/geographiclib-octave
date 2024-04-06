function M = rotm(geog)
% input = phi, lam = geographical coords of viewing direction (deg)
% output = rotation matrix to local cartesian with spherical north

% M = [east; north; up]
  [phi, lam] = deal(geog(1), geog(2));
  M = [-sind(lam), cosd(lam), 0; ...
       -cosd(lam)*sind(phi), -sind(lam)*sind(phi), cosd(phi); ...
       cosd(lam)*cosd(phi), sind(lam)*cosd(phi), sind(phi)];
end
