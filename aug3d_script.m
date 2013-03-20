ACCU = 1;

kind = 0; kindstr = 'helm';

if(1)
  Wall = 2.^[-7:-1];  %Wall = 2.^[-1];
  lowall = cell(length(Wall), 2);
  for g=1:length(Wall)
    W = Wall(g);
    fprintf(1, '%d\n', W);
    tic; tmp = aug3d_lowdata(ACCU,W); toc;
    lowall{g,1} = W;
    lowall{g,2} = tmp;
  end
  
  binstr = sprintf('%s3d_ld_%d.bin',kindstr,ACCU);
  fid = fopen(binstr, 'w');
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
  serialize(fid, lowall, string);
  fclose(fid);
  
  if(0)
    fid = fopen(binstr,'r');
    lowtst = deserialize(fid, string);
    fclose(fid);
  end
end

if(1)
  NPQ = 4; %number of wedges per quadrant
  L = 4096;
  
  Wall = 2.^[0:3];  %Wall = 2.^[0:4];
  hghall = cell(length(Wall), 2);
  res = [];
  for g=1:length(Wall)    %figure;
    W = Wall(g);
    fprintf(1, '%d\n', W);
    respre = res;
    res = aug3d_hghdata(ACCU,NPQ,W,L,respre);    %res = aug3d_hghold(ACCU,NPQ,W,L,respre);
    hghall{g,1} = W;
    hghall{g,2} = res;
  end
  
  binstr = sprintf('%s3d_hd_%d_%d.bin',kindstr,ACCU,NPQ);
  fid = fopen(binstr, 'w');
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
  serialize(fid, hghall, string);
  fclose(fid);
  
  if(0)
    fid = fopen(binstr,'r');
    hghtst = deserialize(fid, string);
    fclose(fid);
  end
end


