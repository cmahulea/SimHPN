%    This is part of SimHPN Toolbox, for Matlab 2010b or newer.
%
%    Copyright (C) 2016 SimHPN developing team. For people, details and citing
%    information, please see: http://webdiis.unizar.es/GISED/?q=tool/simhpn.
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% compute the minimal firing count vector for full rank, or strongly
% connected and consistent structurally persistent PN, in other stuctures
% the minimal firing count vector maybe not unique.

function mfv = SimHPN_control_minFv(C, m0, mf)

[PlaceNum, TransNum] = size(C);


for i = 1 : TransNum
    co(1,i) = 1;
end

lb = zeros(TransNum, 1);
ub = [];

b = mf - m0;

for i = 1 : PlaceNum
   ctype(i,1) = 'S';
end

for i = 1 : TransNum
    vartype(i,1) = 'C';
end

param.msglev = 0;
[xmin, fmin, status, extra] = glpk (co, C, b, lb, ub, ctype, vartype, 1, param);

mfv = xmin;

