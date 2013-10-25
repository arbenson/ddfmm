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
function [points,coords]=new_readwrl(fname, datadir)
%read and normalize
  fprintf('Reading data from %s/%s\n', datadir, fname);
  full_fname = sprintf('%s/%s', datadir, fname);
  fid = fopen(full_fname, 'r');
  assert(fid ~= -1, 'Unsuccessful file open');
  
  tmp = fscanf(fid, '%c', 6);
  points = fscanf(fid, '%f');
  tmp = fscanf(fid, '%c', 6);
  coords = fscanf(fid, '%f');
  
  np = numel(points)/3;
  points = reshape(points, 3, np);
  
  if(coords(4)==-1)
    tmp = 3;
  elseif(coords(5)==-1)
    tmp = 4;
  else
    error('wrong');
  end
  nc = numel(coords)/(tmp+1);
  coords = reshape(coords, tmp+1, nc);
  coords = coords(1:end-1,:);
  coords = coords+1;
  
  if(tmp==4)
    ta = coords(1:3,:);
    tb = coords([4 1 3], :);
    coords = zeros(3,2*size(ta,2));
    coords(:,1:2:end) = ta;
    coords(:,2:2:end) = tb;
  end
  
  %-1 to 1
  ll = min(points, [], 2);
  ur = max(points, [], 2);
  
  md = (ll+ur)/2;
  rr = max(ur-ll)/2;
  
  points = (points-md*ones(1,size(points,2))) / rr;
  %points([2,3],:) = [-points(3,:); points(2,:)];
  
