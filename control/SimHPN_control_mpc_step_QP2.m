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

% the version without the upbound constraints for firing count vector


function [mArray, wArray] = SimHPN_control_mpc_step_QP2( Pre, Post, m0, lambda, delta, mf, wf, Z, Q, R, N )

[placeNum, transNum] = size(Pre);

% construct matrix G
DT = [];
TA = [];

%Z is equal the matrix P for unconstrained LRQ problems
%Z = SimHPN_control_dare_P( Pre, Post, delta, Q, R)

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
% X = [m(k)-mf, m(k+1)-mf,..., m(k+N)-mf, w(k)-wf, w(k+1)-wf,..., w(k+N-1)-wf]'
% k = 0, 1, 2...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VarNum = placeNum * (N+1) + transNum * N;

% constrains matrix for: (m(k+1) - mf) - (m(k) - mf) - C * delta * (w(k) - w_f) = C * delta * wf
CST_1 = [];
CST_2 = [];
CST_1_tmp = [];

for (i = 1 : N)
    CST_1_tmp = blkdiag(CST_1_tmp, eye(placeNum));
    CST_2 = blkdiag(CST_2, delta * C);
end 
%CST_1 = [-1 * CST_1_tmp, zeros(size(CST_1_tmp), placeNum)] + [zeros(size(CST_1_tmp), placeNum), CST_1_tmp];
CST_1 = [-1 * CST_1_tmp, zeros(size(CST_1_tmp,1), placeNum)] + [zeros(size(CST_1_tmp,1), placeNum), CST_1_tmp];
Aeq_1 = [CST_1, -1 * CST_2];
Beq_1 = [];
for i = 1 : N
    Beq_1 = [Beq_1; C * delta * wf];
end

% constrains matrix for: m(k)-mf = m0 - mf
Aeq_2 = [eye(placeNum), zeros(placeNum,N*placeNum + N*transNum)];
Beq_2 = m0 - mf;

relax = 1e-10;

% constrains matrix for: G * [w(k + j) - wf ; m(k + j) - mf] <= G * [-wf; -mf]
CST_3 = [];
for (i = 1 : N)
    temp = zeros(size(G,1), (N + 1)*placeNum + N * transNum);
    
    temp(:, (N+1)*placeNum + (i-1) * transNum + 1 : (N+1)*placeNum + i * transNum) = DT; %G(:, 1:transNum);
    temp(:, (i-1) * placeNum + 1 : (i * placeNum)) = -TA; %G(:, (transNum + 1 : transNum + placeNum));
    
    CST_3 = [CST_3; temp];
end
A_1 = CST_3;
B_1 = [];
for i = 1:N
    B_1 = [B_1; G * [-wf; -mf]];
end

B1 = B_1 + relax * size(B_1,1);


% Low bound 
startIdx = 0;
for i = 1 : N +1
    for j = 1 : placeNum        
        lb(startIdx + placeNum * (i-1) + j, 1) = -mf(j);
    end
end

startIdx = placeNum * (N+1);
for i = 1 : N
    for j = 1 : transNum
        lb(startIdx + transNum * (i-1) + j, 1) = -wf(j)  - relax;
    end
end

% Low bound
ub = Inf * ones(VarNum,1);

H = [];
% objective function
for i = 1 : N
    H = blkdiag(H, Q);
end
H = blkdiag(H, Z);
for i = 1 : N
    H = blkdiag(H, R);
end

f = zeros(VarNum, 1);

Aeq = [Aeq_1; Aeq_2];
Beq = [Beq_1; Beq_2];

A = A_1;
B = B_1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%using quadprog%%%%%%%%%%%%%
X0 = rand((N+1)*placeNum + N * transNum, 1);
X0(1:placeNum,1) = m0 - mf;
options = optimset('Display', 'OFF');
%options=optimset('TolCon', 1e-004);
[X,FVAL,EXITFLAG] = quadprog(H,f,A,B,Aeq,Beq,lb,ub, X0, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mArray = zeros(placeNum, N);
for i = 1 : N
    mArray(:, i)= X(i*placeNum +1 : i*placeNum + placeNum, 1) + mf;
end

wArray = zeros(transNum, N);
startIdx = placeNum * (N+1);
for i = 1 : N
    wArray(:, i) = X(startIdx + (i-1)*transNum +1 : startIdx + (i-1)*transNum + transNum, 1) + wf;
end


if EXITFLAG == -2
    fprintf('***********************************************************\n');
    %m0 = m0
    %mf = mf
    %SimHPN_control_minFv(Post-Pre, m0, mf)     
end

