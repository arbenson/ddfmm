function M = kernel3d(t,s,K)
  
  m = size(t,2);
  n = size(s,2);
  
  tx = t(1,:); ty = t(2,:); tz = t(3,:);
  sx = s(1,:); sy = s(2,:); sz = s(3,:);
  
  x = tx'*ones(1,n) - ones(m,1)*sx;
  y = ty'*ones(1,n) - ones(m,1)*sy;
  z = tz'*ones(1,n) - ones(m,1)*sz;
  
  r = sqrt(x.*x + y.*y + z.*z);
  M = exp(i*K*r) ./ (r);
  
