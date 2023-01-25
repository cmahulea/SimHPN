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

function A = SimHPN_computeA(Pre,Post,lambda)
% for a JF contPN computes the dynamical matrix A

A = zeros(size(Post,1),size(Post,1));
for i = 1 : size(Post,1)
    tout = find(Pre(i,:));
    if ~isempty(tout)
        A(i,i) = -lambda(tout);
    end
    tin = find(Post(i,:));
    for j = 1 : length(tin)
        p = find(Pre(:,tin(j)));
        A(i,p) = (lambda(tin(j)) * Post(i,tin(j))) / Pre(p,tin(j));
    end
end
