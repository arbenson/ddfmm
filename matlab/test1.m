close all;
if(1)
  hold on;
  
  [points, coords] = new_readwrl(fname);
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
  
  curerr = 1e20;
  for it=1:20
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
        if(curwgts(ccpu) < avewgt*1.05)
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
    
    if(0)
      ttl = 0;
      for g=1:NCPU
        cc = cs(:,g);      gud = find(is==g);
        txs = xs(:,gud);      tws = ws(:,gud);
        tmp = ((txs(1,:)-cc(1)).^2 + (txs(2,:)-cc(2)).^2 + (txs(3,:)-cc(3)).^2) .* tws;
        ttl = ttl + sum(tmp);
      end
      fprintf(1, '%d %d\n', it, ttl);
      if(ttl>=0.9999*curerr); break; end;
      curerr = ttl;
    end
     
    colors = 'bgrcmykbgrcmykbgrcmykbgrcmyk'
    for g=1:NCPU
      gud = find(is==g);
      tws = ws(:,gud);
      fprintf(1, '%d %d\n', g, sum(tws));
      plot3(xs(1,gud), xs(2,gud), xs(3,gud), strcat(colors(g),'+'));      axis equal; 
    end
    pause;
  end
  
  
  fprintf(1, 'patch weights\n');
  colors = 'bgrcmykbgrcmykbgrcmykbgrcmyk'
  for g=1:NCPU
    gud = find(is==g);
    tws = ws(:,gud);
    fprintf(1, '%d %d\n', g, sum(tws));
    plot3(xs(1,gud), xs(2,gud), xs(3,gud), strcat(colors(g),'+'));      axis equal; 
  end
end


