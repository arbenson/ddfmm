% Distributed Directional Fast Multipole Method
%   Copyright (C) 2013 Austin Benson, Lexing Ying, and Jack Poulson
%
% This file is part of DDFMM.
%
%    DDFMM is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    DDFMM is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with DDFMM.  If not, see <http://www.gnu.org/licenses/>.
function [prtn, geom] = new_data(fname, datadir, K, NPW, NCPU, NC)

  outdir = sprintf('%s/%s_%d_%d_%d',datadir,fname,K,NPW,NCPU);
  mkdir(outdir);

  [points, coords] = new_readwrl(fname, datadir);
  points = points * K/2 * 0.875; %LEXING: SCALING
  
  %1. get the partition
  numxs = size(coords,2);
  xs = zeros(3,numxs);
  ws = zeros(1,numxs);
  for k=1:numxs
    idx = coords(:,k);
    p1 = points(:,idx(1));
    p2 = points(:,idx(2));
    p3 = points(:,idx(3));  
    xs(:,k) = (p1+p2+p3)/3;
    a = norm(p1-p2);  b = norm(p2-p3);  c = norm(p3-p1);  s = (a+b+c)/2;
    area = sqrt( s*(s-a)*(s-b)*(s-c) );
    npk = ceil(area*NPW*NPW); %total number of points
    ws(k) = npk;
  end
  
  is = zeros(1,numxs);
  
  tmp = floor(rand(1,NCPU) * (numxs-NCPU)) + [1:NCPU];
  cs = xs(:,tmp);

  % Lloyd's Algorithm for point assignment
  for it=1:4*NCPU
    ds = zeros(numxs, NCPU);
    for g=1:NCPU
      cc = cs(:,g);
      ds(:,g) = sqrt((xs(1,:)-cc(1)).^2 + (xs(2,:)-cc(2)).^2 + (xs(3,:)-cc(3)).^2);
    end
    
    is = -ones(1,numxs);
    ttlwgt = sum(ws);
    avewgt = ttlwgt / NCPU;
    curwgts = zeros(1,NCPU);
    [tmp,idxs] = sort(ws);
    for k=numxs:-1:1
      cidx = idxs(k);
      cwgt = ws(cidx);
      dist = ds(cidx,:);
      [dist,order] = sort(dist);
      for g=1:NCPU
        ccpu = order(g);
        if(curwgts(ccpu) <= avewgt*1.05)
          is(cidx) = ccpu;
          curwgts(ccpu) = curwgts(ccpu) + cwgt;
          break;
        end
      end
    end
    if(numel(find(is==-1))>0) error('wrong'); end;
    
    for g=1:NCPU
      gud = find(is==g);      %fprintf(1, '%d\n', numel(gud));
      txs = xs(:,gud);
      tws = ws(:,gud);
      cs(1,g) = sum(txs(1,:).*tws) / sum(tws);
      cs(2,g) = sum(txs(2,:).*tws) / sum(tws);
      cs(3,g) = sum(txs(3,:).*tws) / sum(tws);
    end
    
    colors = 'bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk';
    clf; hold on;
    %for g=1:NCPU
    %  gud = find(is==g);
    %  tws = ws(:,gud);
      %fprintf(1, '%d %d\n', g, sum(tws));
      %view(3);      plot3(xs(1,gud), xs(2,gud), xs(3,gud), strcat(colors(min(g, length(colors))),'+'));      axis equal; 
    %end
    pause(0.2);

    fprintf('iteration number %d\n', it);
  end

  fprintf('end first big loop\n');
  
  % get the prtn
  prtn = cell(NCPU+1,1);
  prtn{1} = 0;
  ttl = 0;
  for g=1:NCPU
    gud = find(is==g);
    tws = ws(:,gud);
    prtn{g+1} = ttl+sum(tws);
    ttl = ttl+sum(tws);
  end
  
  % get the geomptrn
  coef = zeros(NC,NC,NC,NCPU);
  for k=1:numxs
    tmp = floor((xs(:,k) + K/2) / (K/NC)) + 1;
    coef(tmp(1),tmp(2),tmp(3),is(k)) = coef(tmp(1),tmp(2),tmp(3),is(k)) + ws(k);
  end

  coef = reshape(coef, [NC*NC*NC, NCPU]);
  aws = sum(coef,2);
  geom = -ones(1,NC*NC*NC);
  %avewgt = NC*NC*NC / NCPU;
  avewgt = ceil(sum(aws>0)) / NCPU;
  curwgts = zeros(1,NCPU);
  [~, idxs] = sort(aws);
  for k=numel(geom):-1:1
    cidx = idxs(k);
    if(aws(cidx)==0)
      break;
    end
    [~, order] = sort(coef(cidx,:));
    for g=NCPU:-1:1
      ccpu = order(g);
      if(curwgts(ccpu) <= avewgt*1)
        geom(cidx) = ccpu;
        curwgts(ccpu) = curwgts(ccpu) + 1;
        break;
      end
    end
  end
  geom = reshape(geom, [NC,NC,NC]);
  empty = find(geom==-1);
  CPU_ind = 1
  for k = 1:length(empty)
    geom(empty(k)) = CPU_ind;
    CPU_ind = CPU_ind + 1;
    if CPU_ind == NCPU + 1
      CPU_ind = 1;
    end
  end
  geom

  fprintf(1, 'cell weights\n');
  for g=1:NCPU
    fprintf(1, '%d %d\n', g, curwgts(g));
  end
  
  geom = geom - 1; % difference between matlab and c
  binstr = sprintf('%s/geomprtn',outdir);
  string = {'IntNumTns'};
  fid = fopen(binstr, 'w');
  serialize(fid, geom, string);
  fclose(fid);
  
  ttl = sum(ws);
  chk = floor(rand(20,1) * (ttl-20)) + [1:20]';
  binstr = sprintf('%s/chk',outdir);
  string = {'IntNumVec'};
  fid = fopen(binstr, 'w');
  serialize(fid, chk, string);
  fclose(fid);
  
  ttl = 0;
  for g=1:NCPU
    gud = find(is==g);
    tws = ws(:,gud);
    
    pts = zeros(3,30*1000*1000);
    cnt = 0;
    
    %get all the points;
    lib = rand(2,5000);
    good = find(sum(lib,1)<1);
    lib = lib(:,good);
    lib = [lib; 1-sum(lib,1)];
    for k=1:numel(gud)
      idx = coords(:,gud(k));
      p1 = points(:,idx(1));
      p2 = points(:,idx(2));
      p3 = points(:,idx(3));  
      a = norm(p1-p2);  b = norm(p2-p3);  c = norm(p3-p1);  s = (a+b+c)/2;
      area = sqrt( s*(s-a)*(s-b)*(s-c) );
      npk = ceil(area*NPW*NPW);
      tmp = ceil(rand(1,npk) * size(lib,2));
      tmp = lib(:,tmp); %baricenteric
      rsp = p1*tmp(1,:) + p2*tmp(2,:) + p3*tmp(3,:);
      pts(:, cnt+[1:size(rsp,2)]) = rsp;
      cnt = cnt + size(rsp,2);
    end
    
    pts = pts(:,1:cnt);
    den = randn(1,size(pts,2)) + i*rand(1,size(pts,2));
    
    ptsmap = cell(size(pts,2),2);
    for k=1:size(ptsmap,1)
      ptsmap{k,1} = (k-1) + ttl;
      ptsmap{k,2} = pts(:,k);
    end
    
    denmap = cell(size(pts,2),2);
    for k=1:size(denmap,1)
      denmap{k,1} = (k-1) + ttl;
      denmap{k,2} = den(:,k);
    end
    
    %dump data
    binstr = sprintf('%s/pos_%d_%d',outdir,g-1,NCPU);
    string = {'tuple'...
              {'map'...
               {'int'}...
               {'Point3'}...
              }...
              {'vector'...
               {'int'}...
              }...
             };
    fid = fopen(binstr,'w');
    serialize(fid, {ptsmap, prtn}, string);
    fclose(fid);
    
    binstr = sprintf('%s/den_%d_%d',outdir,g-1,NCPU);
    string = {'tuple'...
              {'map'...
               {'int'}...
               {'cpx'}...
              }...
              {'vector'...
               {'int'}...
              }...
             };
    fid = fopen(binstr,'w');
    serialize(fid, {denmap, prtn}, string);
    fclose(fid);
    
    ttl = ttl+sum(tws);
  end
