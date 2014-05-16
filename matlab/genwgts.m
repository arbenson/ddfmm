if(1)
  %only 5 deg now
  n = 5;
  rt = sqrt(15);
  ks = [1/3 (6+rt)/21 (9-2*rt)/21 (6+rt)/21 (6-rt)/21 (9+2*rt)/21 (6-rt)/21 ]';
  es = [1/3 (6+rt)/21 (6+rt)/21 (9-2*rt)/21 (6-rt)/21 (6-rt)/21 (9+2*rt)/21 ]';
  gs = 1-ks-es;
  t1 = (155+rt)/2400;  t2 = (155-rt)/2400;
  ws = [9/80 t1 t1 t1 t2 t2 t2]';
  
  mat = [gs ks es ws*2];
  
  gauwgts = cell(1,2);
  gauwgts{1,1} = n;
  gauwgts{1,2} = mat;
  
  binstr = sprintf('gauwgts.bin');
  fid = fopen(binstr,'w');
  string = {'map', ...
            {'int'}, ...
            {'DblNumMat'}, ...
           };
  serialize(fid, gauwgts, string);
  fclose(fid);
  
  fid = fopen(binstr,'r');
  new = deserialize(fid,string);
  fclose(fid);
end

if(1)
  ns = [4:6];
  sigwgts = cell(length(ns),2);
  for ni=1:numel(ns)
    n = ns(ni);
    if(1)
      beta = 0.5./sqrt(1-(2*(1:n-1)).^(-2));
      T = diag(beta,1) + diag(beta,-1);
      [V,D] = eig(T);
      x = diag(D);
      [x,i] = sort(x);
      w = 2*V(1,i).^2; w = w';
      x = (x+1)/2;
      w = w/2;
    end
    [us,vs] = ndgrid(x);
    ts = w*w';
    as = 1-us;
    bs = us.*(1-vs);
    cs = us.*vs;
    ws = ts.*us;
    mat = [as(:) bs(:) cs(:) 2*ws(:)];
    
    sigwgts{ni,1} = n;
    sigwgts{ni,2} = mat;
  end
  binstr = sprintf('sigwgts.bin');
  fid = fopen(binstr,'w');
  string = {'map', ...
            {'int'}, ...
            {'DblNumMat'}, ...
           };
  serialize(fid, sigwgts, string);
  fclose(fid);
  
  fid = fopen(binstr,'r');
  new = deserialize(fid,string);
  fclose(fid);
end