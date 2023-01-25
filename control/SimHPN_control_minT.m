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

function [s, m, u, w] = SimHPN_control_minT( Pre, Post, m0, lambda, delta, mf, N )

warning off
qdis = 1e-6;
s = 0;
m = [];
u = [];
w = [];

[m, w] = SimHPN_control_minT_LP(Pre, Post, m0, lambda, delta, mf, N);
dif = m(:, N) - mf;
difQ = dif' * dif;  
if difQ > qdis
    errordlg('Please increase the intial guess of steps number OR increase sampling time.','SimHPN: Data error');
    s = 0;
    m = [];
    w = [];
    return;
end


while(1)
    N = N - 1;
   
    [mTmp, cfTmp] = SimHPN_control_minT_LP( Pre, Post, m0, lambda, delta, mf, N );
    dif = mTmp(:, N) - mf;
    difQ = dif' * dif;  
    if difQ <= qdis
        m = mTmp;
        w = cfTmp;
        fprintf('-');
    else
        N = N + 1;
        break;
    end   
end
fprintf('\n');
s = N;
m = [m0, m];

warning on









