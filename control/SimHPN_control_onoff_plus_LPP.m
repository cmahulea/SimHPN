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

function [m, u, w] = SimHPN_control_onoff_plus_LPP( Pre, Post, m, lambda, delta, fv_remain )

%%%%%%%%%%%%%%%%%%%Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = [w_{k}, m_{k+1}]'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% incidence matrix
C = Post - Pre;

[placeNum, transNum] = size(Pre);

% construct matrix G
DT = [];
TA = [];

for (i = 1:placeNum)
    for (j = 1:transNum)
        if (Pre(i, j) > 0)
            tmpD = zeros(1, transNum);
            tmpD(j) = 1;
            tmpTA = zeros(1, placeNum);
            tmpTA(i) = lambda(j)/Pre(i, j);
            DT = [DT;tmpD];
            TA = [TA;tmpTA];
        end
    end
end

G = [DT -TA];

% propositional firing for transitions in confilic
cm = SimHPN_control_cflBalance(Pre);
ctmax = SimHPN_control_conMatx(cm, fv_remain);

% constraints m_k+1 - C * w_k * theta = m_k
A1 = [- C * delta, eye(placeNum)];
b1 = m;
    
% w_k <= remain fv
A2 = [eye(transNum), zeros(transNum, placeNum)];
b2 = fv_remain / delta;

% flow constrains : DT * wk <= TA * mk
A3 = [DT, zeros(size(DT,1), placeNum)];    
b3 = TA * m;

% proportionally fired
if size(ctmax,1) ~= 0 
    A4 = [ctmax, zeros(size(ctmax,1), placeNum)];                
else
    A4 = [];   
end
b4 = zeros(size(A4,1), 1);


A = [A1;A2;A3;A4];
b = [b1;b2;b3;b4];

lb = zeros(placeNum + transNum, 1);
ub = [];

c = [ones(1, transNum), zeros(1, placeNum)];

% sense of each constraint
for i = 1 : placeNum
    ctype1(i, 1) = 'S';
end

for i = 1 : transNum
    ctype2(i, 1) = 'U';
end

for i = 1 : size(TA,1)
    ctype3(i, 1) = 'U';
end

for i = 1 : size(A4, 1)
    ctype4(i,1) = 'S';
end

if size(ctmax,1) ~= 0
    ctype = [ctype1; ctype2; ctype3; ctype4];
else
    ctype = [ctype1; ctype2; ctype3];
end

% type of variables
for i = 1 : placeNum + transNum
    vartype(i,1) = 'C';
end

% sense, maximize
sense = -1; 
param.msglev = 1;

% invoke glpk
[xopt, fmin, status, extra] = glpk (c, A, b, lb, ub, ctype, vartype, sense, param);

m = xopt(transNum+1: transNum+placeNum, 1);
w = xopt(1:transNum, 1);
u = w; %%% not correct;

    
  




