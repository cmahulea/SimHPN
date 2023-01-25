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


function [measuredplaces,measuredcost]=SimHPN_Generate(G,A,F,cost)
% function that makes the tree search looking for the optimal
% cost of the observability

measuredplaces = [];
measuredcost = sum(cost) + 1;


h = 0;
node.v = ones(1,length(G));
node.h = h;
I = [];
nodes{1} = node;

while ~isempty(nodes)
    tempNode = nodes{1};
    current = tempNode.v;
    h = tempNode.h;
    nodes(1)=[];
    cont = 0;
    for i = 1 : size(I,1)
        if (min(current - I(i,:)) >= 0)
            cont = 1;
            break;
        end
    end
    if (cont == 1)
        continue;
    end

    places = [];
    for i = 1 : length(current)
        temp = G{i};
        places(i) =  temp(current(i));
    end
    places = setdiff(places,0);
    if SimHPN_isSolution(union(places,F),A)
        if sum(cost(places)) < measuredcost
            measuredcost = sum(cost(places));
            measuredplaces = places;
        end
    else
        I(size(I,1)+1,:) = current;
        continue;
    end
    next = SimHPN_Children(tempNode,G);
    for i = 1 : length(next)
        nodes{length(nodes)+1} = next{i};
    end
end


