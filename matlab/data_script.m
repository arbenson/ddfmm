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
