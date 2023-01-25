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

function marking = SimHPN_control_mrkBordLpp(PiRow, Pre, Post, m0, mf)

[placeNum, transNum] = size(Pre);

% objective function coefficients
co = zeros(1, placeNum + 1);
co(1, placeNum + 1) = 1;

for i = 1 : (placeNum + 1)
    if (i <= placeNum )
        lb(i, 1) = min(mf(i), m0(i));
    else
        lb(i, 1) = 0;
    end
end

for i = 1 : (placeNum + 1)
    if (i <= placeNum )
        ub(i, 1) = max(mf(i), m0(i)); %mf(i);
    else
        ub(i, 1) = 1;
    end
end

for i = 1 : (placeNum + 1)
    if i == 1
        b(i, 1) = 0;
    else
        b(i, 1) = m0(i - 1);
    end
end


%PiRow
%m0-mf
% constraints coefficients.
A = [PiRow, 0;eye(placeNum), m0-mf];

% the sense of each constraint in the constraint matrix
for i = 1 : placeNum + 1
   ctype(i,1) = 'S';
end

% types of the variables
for i = 1 : placeNum + 1
    vartype(i,1) = 'C';
end

param.msglev = 0;
[xopt, fopt, status, extra] = glpk (co, A, b, lb, ub, ctype, vartype, 1, param);

if status == 5
    marking = xopt(1:placeNum, 1);
else
    marking = 0;
end

