function res = aug3d_hghdata(ACCU,NPQ,W,L,respre)
  
  TPK = 2*pi;
  
  D = W^2+W;  DP = (2*W)^2+2*W;  %D = 2*W*W;  DP = 2*(2*W)*(2*W);
  C = NPQ*W;
  B = C/2;
  
  %get all the centers
  tmp = [0:W:5*D];  %tmp = [-5*D:W:5*D];
  [sx,sy,sz] = ndgrid(tmp);
  sx = sx(:)';  sy = sy(:)';  sz = sz(:)';
  sr = sqrt(sx.^2 + sy.^2 + sz.^2);
  sfts = [sx; sy; sz];
  good = find(sr>=D-1e-8 & sr<=5*D);
  sfts = sfts(:,good);
  
  %the directions of all the centers
  ns = sfts;
  tmp = sqrt(sum(ns.^2,1));
  ns = ns./ (ones(3,1)*tmp);
  ns = aug3d_hghdata_dir(ns, C); %ns = new_nearest(ns, C);
  
  %get the center of the wedges
  [ux,uy,uz] = ndgrid([1:2:C], [1:2:C], C);
  tmp = [ux(:)'; uy(:)'; uz(:)'];
  tmp = sort(tmp);
  uns = unique(tmp', 'rows')';
  
  switch ACCU
   case 1
    EPS = 1e-4;      rr = 70;       RT = 4; NT = 2400;
   case 2
    EPS = 1e-6;      rr = 120;      RT = 5; NT = 2400;
   case 3
    EPS = 1e-8;      rr = 160;      RT = 6; NT = 4800;
  end
  
  res = cell(size(uns,2),2);
  
  for g=1:size(uns,2)
    n = uns(:,g);
    fprintf(1, '%d ', n);    fprintf(1, '\n');
    
    %get pspre and ptpre
    if(W==1)
      if(ACCU==1)
        P = 4;
      elseif(ACCU==2)
        P = 6;
      elseif(ACCU==3)
        P = 8;
      else
        error('wrong ACCU');
      end
      [x,y,z] = ndgrid(linspace(-W/4,W/4,P));
      x = x(:)';  y = y(:)';  z = z(:)';      pp = [x; y; z];
      good = find(max(abs(pp))==W/4);
      pg = pp(:,good);
      pspre = pg;
      ptpre = pg;
    else
      tmp = aug3d_hghdata_predir(n);
      found = 0;
      for k=1:size(respre,1)
        if(prod(double(tmp==respre{k,1})))
          found=found+1;
          idx = k;
        end
      end
      if(found~=1) error('wrong'); end;
      pspre = respre{idx,2}{1}; %uep
      ptpre = -respre{idx,2}{1}; %dcp
    end
    
    %------------- get source and target
    [sx,sy,sz] = ndgrid([-W/4 W/4]); %children box is of size W/2
    sx = sx(:)';  sy = sy(:)';  sz = sz(:)';
    sgs = [sx; sy; sz];
    pscol = [];
    for t=1:size(sgs,2)
      tmp = pspre;      tmp(1,:) = tmp(1,:) + sgs(1,t);      tmp(2,:) = tmp(2,:) + sgs(2,t);      tmp(3,:) = tmp(3,:) + sgs(3,t);
      pscol = [pscol tmp];
    end
    pscol = unique(pscol', 'rows')';
    ptcol = [];
    for t=1:size(sgs,2)
      tmp = ptpre;      tmp(1,:) = tmp(1,:) + sgs(1,t);      tmp(2,:) = tmp(2,:) + sgs(2,t);      tmp(3,:) = tmp(3,:) + sgs(3,t);
      ptcol = [ptcol tmp];
    end
    ptcol = unique(ptcol', 'rows')';
    
    %source
    ps = pscol;
    gg = convhulln(ps');
    salw = unique(gg(:)'); %the important ones    %salw = 1:size(pscol,2);
    
    %target
    tol = 1.0 + 1e-8;
    good = find( abs(ns(1,:)-n(1))<=tol & abs(ns(2,:)-n(2))<=tol & abs(ns(3,:)-n(3))<=1e-8);
    sss = sfts(:,good);
    play = [];
    [sx,sy,sz] = ndgrid([-W/2 W/2]);    sx = sx(:)';  sy = sy(:)';  sz = sz(:)';    sgs = [sx;sy;sz];
    for k=1:size(sss,2)
      if(norm(sss(:,k))<DP-1e-8)
        play = [play 1*sgs + sss(:,k)*ones(1,8)];
      else
        play = [play 3*sgs + sss(:,k)*ones(1,8)];  %LEXING IMPORTANT
      end
    end
    if(1)
      gg = convhulln(play');
      pt0 = play(:,gg);
      talw0 = pt0;
      
      dis = sqrt(sum(sss.^2,1));
      [mindis,cls] = min(dis);
      good = find(dis<mindis+1e-8); %LEXING: important
      pt1 = [];
      for k=good
        pt1 = [pt1 ptcol + sss(:,k)*ones(1,size(ptcol,2))]; %the near ones
      end
      pt1 = unique(pt1', 'rows')';
      talw1 = pt1;
      
      tmp = sqrt(sum(play.^2,1));
      tmp = play ./ (ones(3,1)*tmp);
      thes = atan2(tmp(1,:),tmp(3,:));
      phis = atan2(tmp(2,:),tmp(3,:));
      the_rng = [min(thes), max(thes)];      %mid = [the_rng(1)+the_rng(2)]/2;      rad = [the_rng(2)-the_rng(1)]/2;      the_rng = [mid-rad*1.125 mid+rad*1.125];
      phi_rng = [min(phis), max(phis)];      %mid = [phi_rng(1)+phi_rng(2)]/2;      rad = [phi_rng(2)-phi_rng(1)]/2;      phi_rng = [mid-rad*1.125 mid+rad*1.125];
      %check
      a = n(1);      b = n(2);
      if(the_rng(1)>(a-1)/C*pi/4 | the_rng(2)<(a+1)/C*pi/4 | phi_rng(1)>(b-1)/C*pi/4 | phi_rng(2)<(b+1)/C*pi/4)
        error('wrong');
      end
      [aa,bb] = ndgrid(the_rng,phi_rng);
      %fprintf(1,'%d %d %d %d\n', the_rng(1),the_rng(2),phi_rng(1),phi_rng(2));
      rho_rng = [mindis-W/2 6*D];      %rho_rng = [mindis-W/2 L];
      the = rand(1,NT) * (the_rng(2)-the_rng(1)) + the_rng(1);
      phi = rand(1,NT) * (phi_rng(2)-phi_rng(1)) + phi_rng(1);
      rho = rand(1,NT) * (rho_rng(2)-rho_rng(1)) + rho_rng(1);
      %rho1 = abs(randn(1,NT/2))*4*D; rho2 = rand(1,NT/2)*(rho_rng(2)-rho_rng(1)); rho=[rho1 rho2]+rho_rng(1); %LEXING ???????????
      tmp = sqrt(tan(the).^2 + tan(phi).^2 + 1);
      x = tan(the)./tmp .* rho;
      y = tan(phi)./tmp .* rho;
      z = 1./tmp .* rho;
      pt2 = [x; y; z];
      gg = convhulln(pt2');
      talw2 = pt2(:,unique(gg(:)'));
      
      pt = [pt0 pt1 pt2];
      pt = unique(pt', 'rows')';
      talw = [talw0 talw1 talw2];
      talw = unique(talw', 'rows')';
      [c,ia,ib] = intersect(pt',talw', 'rows');
      talw = ia'; %get the right index
    end
    %fprintf(1,'%d %d\n',length(salw),length(talw));
    
    %------------- generate low rank rep
    M = size(pt,2);
    N = size(ps,2);
        
    sz1 = rr*RT;
    sz2 = rr*RT^2;
    if(1)
      nr = round(sz1*sqrt(M/N));
      if(nr>M)
        nr = M;      rs = [1:M];
      else
        rs = sort(ceil(rand(1,nr)*(M-nr))) + [1:nr];
      end
      rs = unique([rs talw]);
      M1 = kernel3d(pt(:,rs), ps, TPK);
      fprintf(1,'%d %d\n',numel(rs),size(ps,2));
      [Q1,R1,E1] = qr(M1,0); clear M1;    %[Q1,R1,E1] = aug3d_qr(M1,EPS); clear M1;
      gd = find(abs(diag(R1))>EPS*abs(R1(1)));
      Q1 = Q1(:,gd);
      R1 = R1(gd,gd);
      idx1 = E1(gd);    psidx = ps(:,idx1);
      M1 = kernel3d(pt,psidx,TPK);
      [Q1,R1] = qr(M1,0); clear M1;
      
      nc = sz1;
      if(nc>N)
        nc = N;      cs = [1:N];
      else
        cs = sort(ceil(rand(1,nc)*(N-nc))) + [1:nc];
      end
      cs = unique([cs salw idx1]); %VERY IMPORTANT
      M2 = kernel3d(pt, ps(:,cs), TPK);
      fprintf(1,'%d %d\n',size(pt,2), numel(cs));
      [Q2,R2,E2] = qr(M2',0); clear M2;    %[Q2,R2,E2] = aug3d_qr(M2',EPS); clear M2;
      gd = find(abs(diag(R2))>EPS*abs(R2(1)));
      Q2 = Q2(:,gd);
      R2 = R2(gd,gd);
      idx2 = E2(gd);  ptidx = pt(:,idx2);
      M2 = kernel3d(ptidx,ps,TPK);
      [Q2,R2] = qr(M2',0);  clear M2;
    end
    
    if(1)
      nc = sz2;
      if(nc>N)
        nc = N;      cs = [1:N];
      else
        cs = sort(ceil(rand(1,nc)*(N-nc))) + [1:nc];
      end
      cs = unique([idx1 cs]);
      nr = sz2;
      if(nr>M)
        nr = M;      rs = [1:M];
      else
        rs = sort(ceil(rand(1,nr)*(M-nr))) + [1:nr];
      end
      rs = unique([idx2 rs]);
    end
    M3 = kernel3d(pt(:,rs), ps(:,cs), TPK);
    fprintf(1,'%d %d\n',numel(rs), numel(cs));
    G = pinv(Q1(rs,:)) * M3 * pinv(Q2(cs,:)');    clear M3;
    E1 = inv(R1);    E2 = G;    E3 = inv(R2');    %fprintf(1, 'cond numbers: %d %d %d\n', cond(R1), cond(R2), cond(inv(R1) * (G * inv(R2'))));
    fprintf(1, 'rr M N idx1 idx2: %d %d %d %d %d\n', rr, M, N, length(idx1), length(idx2));
    
    if(0)
      hold on;      plot(aa,bb);      plot(aa',bb');
      plot(thes, phis, 'r+');
    end
    if(0)
      figure;            clf;                hold on;
      plot3(ps(1,:), ps(2,:), ps(3,:), 'r+');
      plot3(pt(1,:), pt(2,:), pt(3,:), 'b+');      %G = convhulln(pt');      trimesh(G, pt(1,:), pt(2,:), pt(3,:), 2*ones(size(pt(1,:))));
      plot3(psidx(1,:), psidx(2,:), psidx(3,:), 'g+');
      plot3(ptidx(1,:), ptidx(2,:), ptidx(3,:), 'k+');
      %plot3(play(1,:), play(2,:), play(3,:), 'g+');      %
      %G = convhulln(play');      trimesh(G, play(1,:), play(2,:), play(3,:), ones(size(play(1,:))));
      xlabel('x');      ylabel('y');      zlabel('z');
      view(3);    axis equal; axis tight;
      axis([-3*D 3*D -3*D 3*D -3*D 3*D]);
    end
    
    %----------------------------------
    %check
    test = 0;
    if(0)
      Z = 1000;
      psrnd = (rand(3,Z)-1/2)*W;
      pa = psrnd;
      
      the = rand(1,NT) * (the_rng(2)-the_rng(1)) + the_rng(1);
      phi = rand(1,NT) * (phi_rng(2)-phi_rng(1)) + phi_rng(1);
      rho = rand(1,NT) * (rho_rng(2)-rho_rng(1)) + rho_rng(1);
      tmp = sqrt(tan(the).^2 + tan(phi).^2 + 1);
      x = tan(the)./tmp .* rho;
      y = tan(phi)./tmp .* rho;
      z = 1./tmp .* rho;
      pb = [x; y; z];
      tmpd = randn(size(pa,2),1);
      t1 = kernel3d(pb, pa, 2*pi) * tmpd;
      FWD = kernel3d(ptidx, pa, 2*pi);
      t2 = kernel3d(pb, psidx, 2*pi) * (E1 * (E2 * (E3 * (FWD * tmpd))));
      fprintf(1, 'Error %d %d\n', norm(t1-t2)/norm(t1), max(abs(t1-t2))/max(abs(t1)));
    end
    
    if(0)
      Z = 100;
      psrnd = (rand(3,Z)-1/2)*W;
      pa = psrnd;
      tmp = sort(abs(ns));
      good = find( abs(tmp(1,:)-n(1))<=tol & abs(tmp(2,:)-n(2))<=tol & abs(tmp(3,:)-n(3))<=1e-8 );
      ssa = sfts(:,good);
      nna = ns(:,good);
      disa = sqrt(sum(ssa.^2,1));
      [tmp, clsa] = min(disa);
      tmpd = randn(size(pa,2),1);
      for t=[clsa] %ceil(rand(1,1)*size(ssa,2))]        %ssa(:,t)
        ptrnd = (rand(3,Z)-1/2)*W;
        %first try
        pb = ptrnd + ssa(:,t)*ones(1,size(ptrnd,2));
        [ttt, pat] = sort(abs(nna(:,t)));
        gsidx = psidx;        gtidx = ptidx;
        gsidx(pat,:) = psidx;        gtidx(pat,:) = ptidx;
        sgn = sign(nna(:,t));        sgn(sgn==0) = 1;
        gsidx = (sgn*ones(1,size(gsidx,2))) .* gsidx;
        gtidx = (sgn*ones(1,size(gtidx,2))) .* gtidx;
        t1 = kernel3d(pb, pa, 2*pi) * tmpd;
        FWD = kernel3d(gtidx, pa, 2*pi);
        t2 = kernel3d(pb, gsidx, 2*pi) * (E1 * (E2 * (E3 * (FWD * tmpd))));
        fprintf(1, '%d %d %d %d\n', ssa(1,t), ssa(2,t), ssa(3,t), norm(t1-t2)/norm(t1));
        %hold on;        plot3(pa(1,:), pa(2,:), pa(3,:), 'r+');        plot3(pb(1,:), pb(2,:), pb(3,:), 'b+');
        %second try
        pb = ptrnd + 10*ssa(:,t)*ones(1,size(ptrnd,2));
        [ttt, pat] = sort(abs(nna(:,t)));
        gsidx = psidx;        gtidx = ptidx;
        gsidx(pat,:) = psidx;        gtidx(pat,:) = ptidx;
        sgn = sign(nna(:,t));        sgn(sgn==0) = 1;
        gsidx = (sgn*ones(1,size(gsidx,2))) .* gsidx;
        gtidx = (sgn*ones(1,size(gtidx,2))) .* gtidx;
        t1 = kernel3d(pb, pa, 2*pi) * tmpd;
        FWD = kernel3d(gtidx, pa, 2*pi);
        t2 = kernel3d(pb, gsidx, 2*pi) * (E1 * (E2 * (E3 * (FWD * tmpd))));
        fprintf(1, '%d %d %d %d\n', ssa(1,t), ssa(2,t), ssa(3,t), norm(t1-t2)/norm(t1));
        %hold on;        plot3(pa(1,:), pa(2,:), pa(3,:), 'r+');        plot3(pb(1,:), pb(2,:), pb(3,:), 'b+');
      end
    end
    
    uep = psidx;
    ucp = ptidx;
    uc2ue = {E1 E2 E3}; %uc2ud = ESS;
    
    res{g,1} = n(:);
    res{g,2} = {uep, ucp, uc2ue};    %res{g,2} = uep;    res{g,3} = ucp;    res{g,4} = uc2ue;
    
    %----------------------------
    %switch
    ps_tmp = ps;    ps = -pt;    pt = -ps_tmp;
    psidx_tmp = psidx;    psidx = -ptidx;    ptidx = -psidx_tmp;
    pscol_tmp = pscol;    pscol = -ptcol;    ptcol = -pscol_tmp;    %ESS = transpose(ESS);
    E1_tmp = E1;    E1 = transpose(E3);    E2 = transpose(E2);    E3 = transpose(E1_tmp);
    
    %----------------------------------
    if(0)
      ptrnd = (rand(3,Z)-1/2)*W;
      pb = ptrnd; %pb = ptcol;
      tmp = sort(abs(ns));
      good = find( abs(tmp(1,:)-n(1))<=tol & abs(tmp(2,:)-n(2))<=tol & abs(tmp(3,:)-n(3))<=1e-8 );
      ssa = sfts(:,good);
      nna = ns(:,good);
      disa = sqrt(sum(ssa.^2,1));
      [tmp, clsa] = min(disa);
      tmpd = rand(size(pa,2),1);
      for t=[clsa] %ceil(rand(1,1)*size(ssa,2))]
        psrnd = (rand(3,Z)-1/2)*W;
        %first try
        pa = psrnd - ssa(:,t)*ones(1,size(psrnd,2));
        [ttt, pat] = sort(abs(nna(:,t)));
        gsidx = psidx;        gtidx = ptidx;
        gsidx(pat,:) = psidx;        gtidx(pat,:) = ptidx;
        sgn = sign(nna(:,t));        sgn(sgn==0) = 1;
        gsidx = (sgn*ones(1,size(gsidx,2))) .* gsidx;
        gtidx = (sgn*ones(1,size(gtidx,2))) .* gtidx;
        t1 = kernel3d(pb, pa, 2*pi) * tmpd;
        FWD = kernel3d(pb, gsidx, 2*pi);        %t2 = FWD * ESS * kernel3d(ptidx, pa, 2*pi) * tmp;
        t2 = FWD * (E1 * (E2 * (E3 * (kernel3d(gtidx, pa, 2*pi) * tmpd))));
        fprintf(1, '%d %d %d %d\n', ssa(1,t), ssa(2,t), ssa(3,t), norm(t1-t2)/norm(t1));
        %seocnd try
        pa = psrnd - 10*ssa(:,t)*ones(1,size(psrnd,2));
        [ttt, pat] = sort(abs(nna(:,t)));
        gsidx = psidx;        gtidx = ptidx;
        gsidx(pat,:) = psidx;        gtidx(pat,:) = ptidx;
        sgn = sign(nna(:,t));        sgn(sgn==0) = 1;
        gsidx = (sgn*ones(1,size(gsidx,2))) .* gsidx;
        gtidx = (sgn*ones(1,size(gtidx,2))) .* gtidx;
        t1 = kernel3d(pb, pa, 2*pi) * tmpd;
        FWD = kernel3d(pb, gsidx, 2*pi);        %t2 = FWD * ESS * kernel3d(ptidx, pa, 2*pi) * tmp;
        t2 = FWD * (E1 * (E2 * (E3 * (kernel3d(gtidx, pa, 2*pi) * tmpd))));
        fprintf(1, '%d %d %d %d\n', ssa(1,t), ssa(2,t), ssa(3,t), norm(t1-t2)/norm(t1));
      end
    end
    
    %---------------------
    dep = psidx;
    dcp = ptidx;
    dc2de = {E1 E2 E3};    %dc2dd = ESS;
    
    %res{g,5} = dep;
    %res{g,6} = dcp;
    %res{g,7} = dc2de;
    
  end
  
  
  