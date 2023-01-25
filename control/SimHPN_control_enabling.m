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

% get the enableing degree of transtion j, under makring m%
% Pi is the place control the flow of tj

function [enab, Pi] = SimHPN_control_enabling(Pre, j, m)

ip = find(Pre(:,j));
mi = m(ip);   
Prei = Pre(:, j);
PreiNonZero = Prei(ip);
wmi = mi./PreiNonZero;
[enab, idx] = min(wmi);
Pi = ip(idx);
