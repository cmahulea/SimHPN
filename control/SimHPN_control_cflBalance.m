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

function cm = SimHPN_control_cflBalance(Pre)

global flag;
global cflMatx;

[placeNum, transNum] = size(Pre);

flag = zeros(1, transNum);

% remove the line in Pre which the place has unique output
Pre2 = [];
for i = 1 : placeNum
    if sum(Pre(i,:) ~= 0) > 1
        Pre2 = [Pre2; Pre(i,:)];
    end
end

% set the non-zero elements to be 1
Pre2 = Pre2 > 0;

if size(Pre2,1) == 0
    cm = [];
    return
end

for i = 1 : transNum
    if sum(Pre2(:, i)) == 0
        flag(i) = 1;
    end
end


cflMatx = zeros(placeNum, transNum);
idx = 0;

for j = 1 : transNum
    if flag(j) == 0
        idx = idx + 1;
        SimHPN_control_cflSetGen(j, Pre2, idx);
    end
end

cm = cflMatx(1:idx, :);
