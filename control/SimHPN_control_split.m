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

% insert intermediate state between m0 and mf until the cost time 
% can not decrease by the minimal amount equal to essilon


%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre, Post: Pre- and post- incidence matrix
% m0:        Initial marking
% lambda:    Transitions firing rates
% t:         The costed time for reaching mf from m0
% mf:        Final marking
% epsilon:   Stopping threadhold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Output parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% time: 1 * (interStatNum+1) matrix. time(s) is the time costed for reaching
%           s_th state from (s-1)_th state
% interStatNum: the number of intermediate states being added
% x:    |T| * (interStatNum + 1) matrix. x(:, s) is the accumulated control 
%       for reaching m(:, s - 1) to m(:, s). Notice, m(:, 0) is m0, m(:,
%       interStatNum + 1) = mf.
% m:    |P| * (interStatNum) matrix. m(:, s) is the s_th intermediate state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function  [time, interStatNum, x, m] = SimHPN_control_split(Pre, Post, lambda, m0, mf, t, w, epsilon)

time = [];
x = [];
m = [];

% add a intermediate state md
[md, t1, t2, w1, w2] = SimHPN_control_interState_pw(Pre, Post, lambda, m0, mf);


if (t1 + t2) - t >= 0
   time = t;
   x = w;
   interStatNum = 0;
   m = [];
   return;
end

if  ((t - (t1 + t2))/t > epsilon) 
    [newT1, statNum1, neww1, m1] = SimHPN_control_split(Pre, Post, lambda, m0, md, t1, w1, epsilon);
    [newT2, statNum2, neww2, m2] = SimHPN_control_split(Pre, Post, lambda, md, mf, t2, w2, epsilon); 
    interStatNum = statNum1 + statNum2 + 1;
    time = [newT1, newT2];
    x = [neww1, neww2];
    m = [m1, md, m2];
else
    m = md;
    interStatNum = 1;
    x = [w1, w2]; 
	time = [t1, t2];
end

