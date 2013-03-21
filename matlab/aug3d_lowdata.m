function res = aug3d_lowdata(ACCU,W)
  
  TPK = 2*pi;
  lambda = 1e-16;

  if(    ACCU==1)
    P = 4;
  elseif(ACCU==2)
    P = 6;
  elseif(ACCU==3)
    P = 8;
  end
  
  res = cell(1,4);
  
  %----------------------------------------
  %upeqnpos
  L = W;
  D = P;
  [x,y,z] = ndgrid(linspace(-L/2,L/2,D));
  x = x(:)';  y = y(:)';  z = z(:)';
  pp = [x; y; z];
  good = find(max(abs(pp))==L/2);
  uep = pp(:,good);
  
  %upchkpos
  L = 3*W;
  D = P;
  [x,y,z] = ndgrid(linspace(-L/2,L/2,D));
  x = x(:)';  y = y(:)';  z = z(:)';
  pp = [x; y; z];
  good = find(max(abs(pp))==L/2);
  ucp = pp(:,good);
  
  %upeqnpos_chd  %uep_chd = uep/2;
  
  %up c2d matrix
  Mce = kernel3d(ucp,uep,TPK);
  
  [ua,sa,va] = svd(Mce,0);
  
  gd = find(diag(sa)>lambda*sa(1));
  ua = ua(:,gd);    sa = sa(gd,gd);    va = va(:,gd);
  uc2ue = {va complex(1./diag(sa)) ua'};  %isreal(uc2ue{2})
  
  %----------------------------------------
  %dneqnpos
  dep = ucp;
  
  %dnchkpos
  dcp = uep;
  
  %dnchkpos_chd  dcp_chd = dcp/2;
  
  %dn c2d matrix
  Mce = kernel3d(dcp,dep,TPK);
  [ub,sb,vb] = svd(Mce,0);  %LEXING: ub = conj(va); sb = sa; vb = conj(u)
  gd = find(diag(sb)>lambda*sb(1));
  ub = ub(:,gd);    sb = sb(gd,gd);    vb = vb(:,gd);
  dc2de = {vb complex(1./diag(sb)) ub'};  %isreal(dc2de{2})
  
  %----------------------------------------
  %up d 2 dn c vector
  ue2dc = cell(7,7,7);
  for a=1:7
    for b=1:7
      for c=1:7
        sa = a-4;        sb = b-4;        sc = c-4;
        if(abs(sa)>1 | abs(sb)>1 | abs(sc)>1)
          sfts = [sa;sb;sc]*W;
          
          step = W/(D-1);
          tmp = [0:D-1, -D:-1] * step;
          [x,y,z] = ndgrid(tmp);
          x = x(:)';  y = y(:)';  z = z(:)';
          pc = [x; y; z];
          pc = sfts * ones(1,size(pc,2)) + pc;
          Moc = kernel3d(pc, [0;0;0], TPK);
          
          tmp = reshape(Moc, [2*D,2*D,2*D]);
          tmp = fftn(tmp);
          ue2dc{a,b,c} = tmp;
          %else
          %ud2dc{a,b,c} = rand(2*D,2*D,2*D) + i*rand(2*D,2*D,2*D);
        end
      end
    end
  end
  
  res{1} = uep;
  res{2} = ucp;
  res{3} = uc2ue;  %res{4} = ud2uc;
  
  %res{4} = dep;
  %res{5} = dcp;
  %res{6} = dc2de;  %res{8} = de2dc;
  
  res{4} = ue2dc;
  

  