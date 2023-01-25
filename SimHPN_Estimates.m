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

function ret = SimHPN_Estimates(F,Pre,Post)
% Compute the maximum AF & JF paths that end in F

ret = [];
for i = 1 : length(F)
    p = F(i);
    ret = [ret p];
    tin = find(Post(p,:));
    if isempty(tin)
        continue;
    end
    if length(tin) == 1
        pin = find(Pre(:,tin));
    else
        pin = [];
    end
    while (length(tin) <= 1 && isempty(intersect(pin,ret)))
        pnew = pin;
        ret = [ret pnew];
        p = pnew;
        tin = find(Post(p,:));
        if length(tin) == 1
            pin = find(Pre(:,tin));
        else
            pin = [];
            break;
        end
    end
end

