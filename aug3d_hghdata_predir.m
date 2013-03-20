function r = aug3d_hghdata_predir(n)
  C = max(abs(n));
  r = floor((n+C-1)/2);
  r = 2*floor(r/2) + 1 - C/2;
  r(3) = n(3)/2;
  
  %  r = 1/2*(n+(-1).^((n-1)/2));
  %  r(3) = n(3)/2;
  %  r = [2*r(1)-1; 2*r(2)+1; 2*r(3)];
