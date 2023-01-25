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

function SimHPN_control_cflSetGen(j, Pre2, idx)

[placeNum, transNum] = size(Pre2);

global flag;
global cflMatx;

cflMatx(idx, j) = 1;
flag(1, j) = 1;

for i = 1 : placeNum
     if Pre2(i, j) > 0
        for k = 1 : transNum
            if k ~= j && Pre2(i, k) > 0 && flag(1, k) == 0
                SimHPN_control_cflSetGen(k, Pre2, idx);
            end
        end
    end
end

