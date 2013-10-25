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
function writedata(fid, D)
  if(    iscell(D))
    tp = 0;
  elseif(isreal(D))
    tp = 1;
  else
    tp = 2;
  end
  
  dims = size(D);
  ndim = length(dims);
  fwrite(fid, tp, 'int');
  fwrite(fid, ndim, 'int');
  fwrite(fid, dims, 'int');
  
  if(    tp==0) %cell
    for k=1:numel(D)
      writedata(fid, D{k});
    end
  elseif(tp==1) %real
    R = real(D);
    fwrite(fid,R,'double');
  elseif(tp==2) %complex
    R = real(D);
    I = imag(D);
    [m,n] = size(D);
    T = zeros(size([R;D]));
    T(1:2:end) = R;
    T(2:2:end) = I;
    fwrite(fid,T,'double');
  end
  
  