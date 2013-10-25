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
function M = kernel3d(t,s,K)
  
  m = size(t,2);
  n = size(s,2);
  
  tx = t(1,:); ty = t(2,:); tz = t(3,:);
  sx = s(1,:); sy = s(2,:); sz = s(3,:);
  
  x = tx'*ones(1,n) - ones(m,1)*sx;
  y = ty'*ones(1,n) - ones(m,1)*sy;
  z = tz'*ones(1,n) - ones(m,1)*sz;
  
  r = sqrt(x.*x + y.*y + z.*z);
  M = exp(i*K*r) ./ (r);
  
