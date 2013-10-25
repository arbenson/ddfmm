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

ACCU = 3;

kind = 0; kindstr = 'helm';
datadir = './..';

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
  
  binstr = sprintf('%s/data/%s3d_ld_%d.bin',datadir,kindstr,ACCU);
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
end

if(1)
  NPQ = 4; %number of wedges per quadrant
  L = 4096;
  
  Wall = 2.^[0:4];  %Wall = 2.^[0:4];
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
  
  binstr = sprintf('%s/data/%s3d_hd_%d_%d.bin',datadir,kindstr,ACCU,NPQ);
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
end
