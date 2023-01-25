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

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mArray:  |P|* N matrix, mArray(:, k) is the state in k_th step
% wArray:  |T|* N matrix, wArray(:, k) is the controlled flow in k_th step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [mArray, wArray] = SimHPN_control_minT_LP( Pre, Post, m0, lambda, delta, mf, N )

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
C = Post - Pre;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X = [m(k)-mf, m(k+1)-mf,..., m(k+N)-mf, w(k), w(k+1)-wf,..., w(k+N-1)-wf, md(k+1), md(k+2),...,md(k+N)]'
% k = 0, 1, 2...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wf = zeros(transNum, 1);
VarNum = placeNum * (N+1) + transNum * N + N * placeNum;

% constrains matrix for: (m(k+1) - mf) - (m(k) - mf) - C * delta * (w(k) - w_f) = C * delta * wf
CST_1 = [];
CST_2 = [];
CST_1_tmp = [];

for (i = 1 : N)
    CST_1_tmp = blkdiag(CST_1_tmp, eye(placeNum));
    CST_2 = blkdiag(CST_2, delta * C);
end 
CST_1 = [-1 * CST_1_tmp, zeros(size(CST_1_tmp), placeNum)] + [zeros(size(CST_1_tmp), placeNum), CST_1_tmp];

Aeq_1 = [CST_1, -1 * CST_2, zeros(size(CST_1,1), N*placeNum)];
Beq_1 = [];
for i = 1 : N
    Beq_1 = [Beq_1; C * delta * wf];
end

% constrains matrix for: m(k)-mf = m0 - mf
Aeq_2 = [eye(placeNum), zeros(placeNum, N*placeNum + N*transNum + N*placeNum)];
Beq_2 = m0 - mf;


% constrains matrix for: G * [w(k + j) - wf ; m(k + j) - mf] <= G * [-wf; -mf]
CST_3 = [];
for (i = 1 : N)
    temp = zeros(size(G,1), (N + 1)*placeNum + N * transNum);
    
    temp(:, (N+1)*placeNum + (i-1) * transNum + 1 : (N+1)*placeNum + i * transNum) = DT; %G(:, 1:transNum);
    temp(:, (i-1) * placeNum + 1 : (i * placeNum)) = -TA; %G(:, (transNum + 1 : transNum + placeNum));
    
    CST_3 = [CST_3; temp];
end
A_1 = [CST_3, zeros(size(CST_3,1), N*placeNum)];
B_1 = [];
for i = 1:N
    B_1 = [B_1; G * [-wf; -mf]];
end


% constrains matrix for: m(k+1) - mf <= md(k+1) and -m(k+1) + mf <= md(k+1)
temp1 = [];
temp2 = [];
for i = 1 : N
    temp1 = blkdiag(temp1, eye(placeNum));
    temp2 = blkdiag(temp2, -1 * eye(placeNum));
end
A_2 = [zeros(size(temp1,1),placeNum), temp1, zeros(size(temp1,1),N*transNum), temp2];
B_2 = zeros(size(A_2,1),1);

A_3 = [zeros(size(temp1,1),placeNum), -1*temp1, zeros(size(temp1,1),N*transNum), temp2];
B_3 = zeros(size(A_3,1),1);

% constrains matrix for: w(k) - mf <= wd(k) and -w(k) + mf <= wd(k)
%temp1 = [];
%temp2 = [];
%for i = 1 : N
%    temp1 = blkdiag(temp1, eye(transNum));
%    temp2 = blkdiag(temp2, -1 * eye(transNum));
%end
%A_4 = [zeros(size(temp1,1),placeNum * (N+1)), temp1, zeros(size(temp1,1),N*placeNum), temp2];
%B_4 = zeros(size(A_4,1),1);

%A_5 = [zeros(size(temp1,1),placeNum * (N+1)), -1*temp1, zeros(size(temp1,1),N*placeNum), temp2];
%B_5 = zeros(size(A_5,1),1);

% w0 + w1 + ... <= fv_remain
%temp = [];
%for i = 1:N
%    temp = [temp, eye(transNum)];
%end
%rows = size(temp,1);
%A_6 = [zeros(rows, (N+1) * placeNum), temp * delta, zeros(rows, N*transNum+N*placeNum)];
%B_6 = fv_remain - N * wf * delta; %????

% w0 > epsilon
%A_7 = [zeros(1, placeNum * (N+1)), -1 * ones(1, transNum), zeros(1, transNum * (N-1)), zeros(1, N*transNum+N*placeNum)];
%B_7 = sum(wf) - 1e-16;

% Low bound 
lb = zeros(VarNum,1);

startIdx = 0;
for i = 1 : N +1
    for j = 1 : placeNum        
        lb(startIdx + placeNum * (i-1) + j, 1) = -mf(j);
    end
end
startIdx = placeNum * (N+1);
for i = 1 : N
    for j = 1 : transNum
        lb(startIdx + transNum * (i-1) + j, 1) = -wf(j);
    end
end

% Low bound
ub = Inf * ones(VarNum,1);

Aeq = [Aeq_1; Aeq_2];
Beq = [Beq_1; Beq_2];

A = [A_1;A_2;A_3];
B = [B_1;B_2;B_3];

c = zeros(1, VarNum);
startIdx = (N+1) * placeNum + N * transNum;
for i = startIdx + (N - 1) * placeNum + 1 : startIdx + N * placeNum
    c(i) = 1;
end

param.msglev = 0;

AA = [Aeq;A];
BB = [Beq;B];
for i = 1 : size(Aeq,1)
    ctype(i,1) = 'S';
end
for i = size(Aeq,1) + 1 : size(Aeq,1)+size(A,1)
    ctype(i,1) = 'U';
end

for i = 1 : VarNum
    vartype(i,1) = 'C';
end

[X, fopt, status, extra] = glpk (c, AA, BB, lb, ub, ctype, vartype, 1, param);
if (status ~= 5)
    fprintf('glpk didnot find optimal solution\n');
    halt
end

mArray = zeros(placeNum, N);
for i = 1 : N
    mArray(:, i)= X(i*placeNum +1 : i*placeNum + placeNum, 1) + mf;
end

wArray = zeros(transNum, N);
startIdx = placeNum * (N+1);
for i = 1 : N
    wArray(:, i) = X(startIdx + (i-1)*transNum +1 : startIdx + (i-1)*transNum + transNum, 1) + wf;
end
