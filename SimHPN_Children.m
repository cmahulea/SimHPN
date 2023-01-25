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


function ret = SimHPN_Children(node,G)
% computes the children of a node


ret = {};
v = node.v;
h = node.h + 1;

while (h <= length(G))
    if (v(h) < length(G{h}))
        new_v = v;
        new_v(h) = new_v(h) + 1;
        new_node.v = new_v;
        new_node.h = h - 1;
        ret{length(ret)+1} = new_node;
    end
    h = h + 1;
end
return

