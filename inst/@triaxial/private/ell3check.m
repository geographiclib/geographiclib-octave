function result = ell3check(ell)
  result = numel(ell) == 3 && size(ell, 1) == 1 && isfinite(ell(1)) && ...
           ell(1) >= ell(2) && ell(2) >= ell(3) && ell(3) > 0;
end
