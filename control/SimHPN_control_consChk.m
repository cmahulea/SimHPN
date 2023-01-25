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

% To check if the net is consistent
% return 1 if the net is consistent, otherwise return 0

function cons = SimHPN_control_consChk(C)


[PlaceNum, TransNum] = size(C);

for i = 1 : TransNum
    co(1,i) = 1;
end

lb = zeros(TransNum, 1);
for i = 1 : TransNum
    lb(i, 1) = 0.01;
end

ub = [];

b = zeros(PlaceNum, 1);

for i = 1 : PlaceNum
   ctype(i,1) = 'S';
end

for i = 1 : TransNum
    vartype(i,1) = 'C';
end

param.msglev = 0;
[xopt, fmin, status, extra] = glpk (co, C, b, lb, ub, ctype, vartype, 1, param);

if status == 5
    cons = 1;
else
    cons = 0;
    %fprintf('The net is not consistent\n');
end

    
