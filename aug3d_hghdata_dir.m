function r = aug3d_hghdata_dir(n, C)
  %map directions to the sphere parameterization used in the paper
  B = C/2;
  t = max(abs(n),[],1);
  n = n./(ones(3,1)*t);
  n = atan(n)/(pi/4);
  r = n*C;
  
  