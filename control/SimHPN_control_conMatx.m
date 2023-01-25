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

% the constraint matrix, fire conflicting transitions proportionally

function CtMatx = SimHPN_control_conMatx(cm, fv)

[rowNum, colNum] = size(cm);

CtMatx = [];
CtNum = 0;

for i = 1 : rowNum
    for j = 1 : colNum - 1
        for k = j + 1 : colNum
            if cm(i, j) > 0 && cm(i, k) > 0 && fv(j) > 0 && fv(k) > 0
                CtNum = CtNum + 1;
                CtMatx(CtNum, j) = 1;
                CtMatx(CtNum, k) = - fv(j)/fv(k);
            end
        end
    end
end

temp = size(CtMatx, 2);
if temp < colNum
    CtMatx = [CtMatx, zeros(size(CtMatx,1), colNum - temp)];
end

