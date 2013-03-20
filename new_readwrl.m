function [points,coords]=new_readwrl(fname)
%read and normalize
  fid = fopen(fname, 'r');
  
  tmp = fscanf(fid, '%c', 6);
  points = fscanf(fid, '%f');
  tmp = fscanf(fid, '%c', 6);
  coords = fscanf(fid, '%f');
  
  np = numel(points)/3;
  points = reshape(points, 3, np);
  
  if(coords(4)==-1)
    tmp = 3;
  elseif(coords(5)==-1)
    tmp = 4;
  else
    error('wrong');
  end
  nc = numel(coords)/(tmp+1);
  coords = reshape(coords, tmp+1, nc);
  coords = coords(1:end-1,:);
  coords = coords+1;
  
  if(tmp==4)
    ta = coords(1:3,:);
    tb = coords([4 1 3], :);
    coords = zeros(3,2*size(ta,2));
    coords(:,1:2:end) = ta;
    coords(:,2:2:end) = tb;
  end
  
  %-1 to 1
  ll = min(points, [], 2);
  ur = max(points, [], 2);
  
  md = (ll+ur)/2;
  rr = max(ur-ll)/2;
  
  points = (points-md*ones(1,size(points,2))) / rr;
  %points([2,3],:) = [-points(3,:); points(2,:)];
  