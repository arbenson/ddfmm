datadir = '../data'

if(0)
  fname = 'F16.wrl';
  K = 64;
  NPW = 20;
  NCPU = 1;
  NC = 8;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end

if(1)
  fname = 'F16.wrl';
  K = 64;
  NPW = 20;
  NCPU = 8;
  NC = 8;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end

if(0)
  fname = 'F16.wrl';
  K = 64;
  NPW = 20;
  NCPU = 16;
  NC = 8;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end

%-----------------------------------
if(0)
  fname = 'SubmarineJ.wrl';
  K = 64;
  NPW = 20;
  NCPU = 4;
  NC = 8;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end

if(0)
  fname = 'SubmarineJ.wrl';
  K = 64;
  NPW = 20;
  NCPU = 4;
  NC = 8;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end

if(0)
  fname = 'SubmarineJ.wrl';
  K = 64;
  NPW = 20;
  NCPU = 1;
  NC = 8;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end

if(0)
  fname = 'SubmarineJ.wrl';
  K = 2;
  NPW = 20;
  NCPU = 1;
  NC = 2;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end

%-----------------------------------
if(0)
  fname = 'sphere.wrl';
  K = 64;
  NPW = 20;
  NCPU = 1;
  NC = 2;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end

if(0)
  fname = 'sphere.wrl';
  K = 64;
  NPW = 20;
  NCPU = 32;
  NC = 16;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end

if(0)
  fname = 'sphere.wrl';
  K = 16;
  NPW = 20;
  NCPU = 1;
  NC = 2;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end

if(0)
  fname = 'sphere.wrl';
  K = 16;
  NPW = 20;
  NCPU = 2;
  NC = 2;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end

if(0)
  fname = 'sphere.wrl';
  K = 16;
  NPW = 20;
  NCPU = 4;
  NC = 4;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end


if(0)
  fname = 'sphere.wrl';
  K = 2;
  NPW = 20;
  NCPU = 1;
  NC = 1;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end

if(0)
  fname = 'sphere.wrl';
  K = 2;
  NPW = 20;
  NCPU = 2;
  NC = 2;
  [prtn, geom] = new_data(fname,datadir,K,NPW,NCPU,NC);
end
