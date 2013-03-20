if(1)
  
end

if(0)
  ACCU = 1;
  NPQ = 4;
  kind = 0; kindstr = 'helm';

  binstr = sprintf('aug3d_ld_%d.bin',ACCU);
  fid = fopen(binstr,'r');
  lowall = readdata(fid);
  fclose(fid);

  binstr = sprintf('aug3d_hd_%d_%d.bin',ACCU,NPQ);
  fid = fopen(binstr,'r');
  hghall = readdata(fid);
  fclose(fid);
end

if(0)
  num = size(lowall,1);
  lownew = cell(num,2);
  for k=1:num
    lownew{k,1} = lowall{k,1};
    lownew{k,2} = lowall(k,2:end);
  end
  
  binstr = sprintf('helm3d_ld_%d.bin',ACCU);
  fid = fopen(binstr,'w');
  string = {'map', ...
            {'double'}, ...
            {'tuple', ...
             {'DblNumMat'}, ...
             {'DblNumMat'}, ...
             {'NumVec', ...
              {'CpxNumMat'} ...
             }, ...
             {'NumTns', ...
              {'CpxNumTns'}...
             }...
            }...
           };
  serialize(fid, lownew, string);
  fclose(fid);
  
  fid = fopen(binstr,'r');
  lowtst = deserialize(fid, string);
  fclose(fid);
end

if(0)
  num = size(hghall,1);
  hghnew = cell(num,2);
  for k=1:num
    hghnew{k,1} = hghall{k,1};
    hghnew{k,2} = cell(size(hghall{k,2},1),2);
    for g=1:size(hghall{k,2},1)
      hghnew{k,2}{g,1} = hghall{k,2}{g,1};
      hghnew{k,2}{g,2} = hghall{k,2}(g,2:end);
    end    
  end
  
  binstr = sprintf('helm3d_hd_%d_%d.bin',ACCU,NPQ);
  fid = fopen(binstr,'w');
  string = {'map' ...
            {'double'} ...
            {'map' ...
             {'Index3'} ...
             {'tuple' ...
              {'DblNumMat'} ...
              {'DblNumMat'} ...
              {'NumVec' ...
               {'CpxNumMat'} ...
              }...
             }...
            }...
           };
  serialize(fid, hghnew, string);
  fclose(fid);
  
  fid = fopen(binstr,'r');
  hghtst = deserialize(fid, string);
  fclose(fid);
end

if(0)
  N = 40;
  [tx,ty,tz] = sphere(N);
  surf(tx,ty,tz);  axis equal;
  shading interp;
  lightangle(-45,30);
  ti = zeros(size(tx));
  ti(:) = 1:numel(ti);
  id = 1:N;
  a1 = ti(id,id);
  a2 = ti(id+1,id);
  a3 = ti(id+1,id+1);
  a4 = ti(id, id+1);
  ll = [a1(:) a2(:) a3(:)]';
  ur = [a1(:) a3(:) a4(:)]';
  points = [tx(:) ty(:) tz(:)]';
  coords = [ll ur];
  trimesh(coords', points(1,:), points(2,:), points(3,:));  axis equal;
  
  coords = coords-1;
  fid = fopen('sphere.wrl','w');
  fprintf(fid, 'points\n');
  for k=1:size(points,2)
    fprintf(fid, '%d %d %d \n', points(1,k), points(2,k), points(3,k));
  end
  fprintf(fid, 'coords\n');
  for k=1:size(coords,2)
    fprintf(fid, '%d %d %d -1 \n', coords(1,k), coords(2,k), coords(3,k));
  end
  fclose(fid);
end


