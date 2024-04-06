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
