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

function [sigmaV, m_steady, flow_steady] = minFv_optCtl(Pre, Post, lambda, m0, tIdx, flowMax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_{ss}, \sigma, m_{ss}, v(.,t1) v(.,t_2)...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[placeNum, transNum] = size(Pre);

% compute the number of variables
varNum = transNum + transNum + placeNum;
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

% constrain 1: C * f_ss = 0
CST_1 = [C, zeros(size(C, 1), varNum - size(C, 2))];

% constrain 2: m = m0 + C * \sigma -> - C * \sigma + m = m0
CST_2 = [zeros(placeNum, transNum), -1 * C, eye(placeNum), zeros(placeNum, vNum)];

% constrain 3: \f_ss[i] = \lambda_i * (m_j / Pre(p_j, t_i) - v(p_j, t_i)
CST_3 = [];
idx = 0;

for ( j = 1 : transNum)
    Prej = Pre(:, j);
    for k = 1 : placeNum
        newRow = [];
        if Prej(k) ~= 0
            idx = idx + 1;
            newRow = zeros(1, varNum);
            newRow(j) = 1;
            newRow(2 * transNum + k) = -lambda(j) / Pre(k, j);
            newRow(2 * transNum + placeNum + idx) = 1;
        end
        CST_3 = [CST_3; newRow];
    end
end

% constrain 4: f_ss[tIdx] = flowMax
CST_4 = zeros(1, varNum);
CST_4(1, tIdx) = 1;


A = [CST_1; CST_2; CST_3; CST_4];


% types of the variables
for i = 1 :  varNum
    vartype(i,1) = 'C';
end

% the sense of each constraint in the constraint matrix
% constrain 1: C * \phi = 0
for i = 1 : size(CST_1, 1)
   b(i, 1) = 0;
   ctype(i,1) = 'S';
end

% constrain 2: m = m0 + C * \sigma -> - C * \sigma + m = m0
for i = size(CST_1, 1) + 1 : size(CST_1, 1) + size(CST_2, 1)
   b(i, 1) = m0(i - size(CST_1, 1));
   ctype(i,1) = 'S';
end

% constrain 3: f_ss = \lambda_i * (m_j / Pre(p_j, t_i) - v(p_j, t_i)
for i = size(CST_1, 1) + size(CST_2, 1) + 1 : size(CST_1, 1) + size(CST_2, 1) + size(CST_3, 1)
   b(i, 1) = 0;
   ctype(i,1) = 'S';
end

% constrain 4:
ctype(size(CST_1, 1) + size(CST_2, 1) + size(CST_3, 1)+1, 1) = 'S';


% constrain 4: f_ss[tIdx] = flowMax
b(size(CST_1, 1) + size(CST_2, 1) + size(CST_3, 1) + 1, 1) = flowMax;

lb = zeros(varNum, 1);
ub = [];


c1 = zeros(1, transNum);
c2 = ones(1, transNum);
%c2 = 1./lambda';
c3 = zeros(1, placeNum);


c = [c1, c2, c3, zeros(1, varNum - 2*transNum - placeNum)];

param.msglev = 0;

[xopt, fopt, status, extra] = glpk (c, A, b, lb, ub, ctype, vartype, 1, param);
if (status ~= 5)
    fprintf('glpk didnot find optimal solution\n');
    halt
end


sigmaV = xopt(transNum + 1 : 2*transNum, 1);
m_steady = xopt(2 * transNum + 1 : 2 * transNum + placeNum, 1);
flow_steady = xopt(1: transNum, 1);

%fopt
%xopt
