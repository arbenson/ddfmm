function D = readdata(fid)
  
  tp = fread(fid, 1, 'int');
  ndim = fread(fid, 1, 'int');
  dims = fread(fid, [1,ndim], 'int');
  
  if(tp==0) %cell
    D = cell(dims);
    for k=1:prod(dims)
      D{k} = readdata(fid);
    end
  elseif(tp==1) %real
    D = fread(fid, prod(dims), 'double');
    D = reshape(D, dims);
  elseif(tp==2) %complex
    T = fread(fid, 2*prod(dims), 'double');
    R = reshape(T(1:2:end), dims);
    I = reshape(T(2:2:end), dims);
    D = R + i*I;
  end
    
  