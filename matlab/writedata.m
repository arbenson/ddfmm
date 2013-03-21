function writedata(fid, D)
  if(    iscell(D))
    tp = 0;
  elseif(isreal(D))
    tp = 1;
  else
    tp = 2;
  end
  
  dims = size(D);
  ndim = length(dims);
  fwrite(fid, tp, 'int');
  fwrite(fid, ndim, 'int');
  fwrite(fid, dims, 'int');
  
  if(    tp==0) %cell
    for k=1:numel(D)
      writedata(fid, D{k});
    end
  elseif(tp==1) %real
    R = real(D);
    fwrite(fid,R,'double');
  elseif(tp==2) %complex
    R = real(D);
    I = imag(D);
    [m,n] = size(D);
    T = zeros(size([R;D]));
    T(1:2:end) = R;
    T(2:2:end) = I;
    fwrite(fid,T,'double');
  end
  
  