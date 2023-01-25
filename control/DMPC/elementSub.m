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

% get the P/T element connected with tj, in the subsystem
function [subNetP, subNetT] = elementSub(Pre, Post, buffer, tj, subNetPold, subNetTold)

subNetP = [];
subNetT = [];

if isempty(find(subNetTold == tj))
    subNetT = [subNetT, tj];    
    
    for i = 1: size(Pre,1)
        if (Pre(i, tj) > 0 || Post(i,tj) > 0) && (isempty(find(buffer == i))) && (isempty(find([subNetPold, subNetP] == i)))
            subNetP = [subNetP, i];
            for j = 1 : size(Pre, 2)
                if Pre(i,j) > 0 || Post(i,j) > 0
                    [pTemp, tTemp] = elementSub(Pre, Post, buffer, j, [subNetP, subNetPold], [subNetT, subNetTold]);
                    subNetP = [subNetP, pTemp];
                    subNetT = [subNetT, tTemp];
                end
            end
        end

    end
end
