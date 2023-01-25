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

% to generate the configuration matrix in a marking m

function Cfg = SimHPN_control_CfgMatrixGen(Pre, Post, m)

[placeNum, transNum] = size(Pre);

for (j = 1:transNum)
    Prej = Pre(:, j);
    nonZero = find(Prej);
    if (size(nonZero, 2) == 0)
        fprintf('trans without input plc\n');
    end
    
    enab = m(nonZero) ./ Prej(nonZero);
    [minV, minIdx] = min(enab);
    plcIdx = nonZero(minIdx);
    Cfg(j, :) = zeros(1, placeNum);
    Cfg(j, plcIdx) = 1 / Prej(plcIdx);
end
