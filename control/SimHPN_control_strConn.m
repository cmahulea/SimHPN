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

% if the net is strongly connected return 1, otherwise return 0

function sc = SimHPN_control_strConn(Pre, Post)

[placeNum, transNum] = size(Pre);

% get the matirx of the conrresponding graph
G = [zeros(placeNum, placeNum), Pre; Post', zeros(transNum, transNum)];

[S, D] = graphconncomp(sparse(G));

%view(biograph(G))

[tempR, tempC] = size(D);
if (S == 1) && (tempC == placeNum + transNum) && isequal(D, ones(1, placeNum + transNum))
    sc = 1;
else
    sc = 0;
    %fprintf('The net is not strongly connected\n');
end


