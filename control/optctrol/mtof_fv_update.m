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

function [w_ss, fv, m_ss, maxStep] = mtof_fv_update(Pre, Post, lambda, delta, m0, tj, optFlow)

%%%%%%%%%%%%%%%%%%%%%%%%%%%Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w_{ss}, fv, m_{ss}, d, v(.,t1) v(.,t_2)...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[placeNum, transNum] = size(Pre);

% compute the number of variables
varNum = transNum + transNum + placeNum + 1;
vNum = 0;

for i = 1 : placeNum
    for j = 1 : transNum
        if Pre(i,j) ~= 0
            varNum = varNum + 1;
            vNum = vNum + 1;
        end
    end
end

C = Post - Pre;
ctype = [];

% constrain 1: C * w_ss = 0
CST_1 = [C, zeros(size(C, 1), varNum - size(C, 2))];
b1 = zeros(placeNum,1);
for i = 1 : placeNum
    ctype = strcat(ctype, 'S');
end

% constrain 2: m = m0 + C * fv -> - C * fv + m = m0
CST_2 = [zeros(placeNum, transNum), -1 * C, eye(placeNum), zeros(placeNum, vNum+1)];
b2 = m0;
for i = 1 : placeNum
    ctype = strcat(ctype, 'S');
end

% constrain 3: w_ss[t] = \lambda_j * (m_ss[p_i] / Pre(p_i, t) - v(p_i, t)
CST_3 = [];
idx = 0;

for j = 1 : transNum
    Prej = Pre(:, j);
    for k = 1 : placeNum
        newRow = [];
        if Prej(k) ~= 0
            idx = idx + 1;
            newRow = zeros(1, varNum);
            newRow(j) = 1;
            newRow(2 * transNum + k) = -lambda(j) / Pre(k, j);
            newRow(2 * transNum + placeNum + 1 + idx) = 1;
        end
        CST_3 = [CST_3; newRow];
    end
end
b3 = zeros(vNum, 1);
for i = 1 : vNum
    ctype = strcat(ctype, 'S');
end

% constraint 4: w_ss[tj] = optFlow
CST_4 = zeros(1, varNum);
CST_4(tj) = 1;
b4 = optFlow;
ctype = strcat(ctype,'S');


% constraint 5: d >= fv(tj) / flow(tj) for any tj

epsilon = 1e-6;
% compute the (uncontrolled) flow vector
flow = zeros(transNum, 1);
for j = 1:transNum       
    % enabling degree
    ip = find(Pre(:,j));        
    mi = m0(ip);   
    Prei = Pre(:, j);
    PreiNonZero = Prei(ip);
    wmi = mi./PreiNonZero;
    enab(j) = min(wmi);        
    flow(j) = lambda(j) * enab(j);
    if flow(j) == 0
        flow(j) = epsilon; % keep the flow to be positive
    end
end   

CST_5 = [];
for j = 1: transNum
    newRow = zeros(1, varNum);
    newRow(1, 2 * transNum + placeNum + 1) = 1;
    newRow(1, transNum + j) = - 1 / flow(j);
    CST_5 = [CST_5; newRow];
end
b5 = zeros(transNum,1);
for i = 1 : transNum
    ctype = strcat(ctype, 'L');
end

A = [CST_1; CST_2; CST_3; CST_4; CST_5];
b = [b1;b2;b3;b4;b5];

% types of the variables
for i = 1 :  varNum
    vartype(i,1) = 'C';
end

lb = zeros(varNum, 1);
ub = [];

% optimizing weight
w1 = 0;
w2 = 10;

c1 = zeros(1,transNum);
c2 = ones(1, transNum) * w1;
%c2(1) = -1;
c3 = zeros(1, placeNum);

c = [c1, c2, c3, w2, zeros(1, vNum)];

param.msglev = 0;

[xopt, fopt, status, extra] = glpk (c, A, b, lb, ub, ctype, vartype, 1, param);
if (status ~= 5)
    fprintf('glpk didnot find optimal solution\n');
    halt
end

w_ss = xopt(1 : transNum, 1);
fv = xopt(transNum + 1 : 2 * transNum, 1);
m_ss = xopt(2*transNum + 1: 2 * transNum + placeNum, 1);
maxStep = xopt(2*transNum+placeNum+1,1)/delta;

%fopt
%xopt
