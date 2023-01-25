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

%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre, Post: Pre- and post- incidence matrix
% m0:        Initial marking
% lambda:    Transitions firing rates
% delta:     Sampling period
% mf:        Final marking
% wf:        Final flow
% fv_remain: The firing count vector
% Z, Q, R:   Weight matices
% N:         Horizon step
% epsilon:   the lower bound of the marking of any place
% zeta:      the lower bound of firing flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mArray:  |P|* N matrix, mArray(:, k) is the state in k_th step
% wArray:  |T|* N matrix, wArray(:, k) is the controlled flow in k_th step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mArray, wArray] = SimHPN_control_mpc_onoff_step_QP2b(Pre, Post, m0, lambda, delta, mf, wf, fv_remain, Z, Q, R, N, epsilon, zeta)

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
% X = [m(k)-mf, m(k+1)-mf,..., m(k+N)-mf, w(k)-wf, w(k+1)-wf,..., w(k+N-1)-wf, relax]'
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
CST_1 = [-1 * CST_1_tmp, zeros(size(CST_1_tmp), placeNum)] + [zeros(size(CST_1_tmp), placeNum), CST_1_tmp];

Aeq_1 = [CST_1, -1 * CST_2, zeros(size(CST_1,1), 1)];
Beq_1 = [];
for i = 1 : N
    Beq_1 = [Beq_1; C * delta * wf];
end

% constrains matrix for: m(k)-mf = m0 - mf
Aeq_2 = [eye(placeNum), zeros(placeNum, N*placeNum + N*transNum + 1) ];
Beq_2 = m0 - mf;

% constrains matrix for: G * [w(k + j) - wf ; m(k + j) - mf] <= G * [-wf; -mf]
CST_3 = [];
for (i = 1 : N)
    temp = zeros(size(G,1), (N + 1)*placeNum + N * transNum);
    
    temp(:, (N+1)*placeNum + (i-1) * transNum + 1 : (N+1)*placeNum + i * transNum) = DT; %G(:, 1:transNum);
    temp(:, (i-1) * placeNum + 1 : (i * placeNum)) = -TA; %G(:, (transNum + 1 : transNum + placeNum));
    
    CST_3 = [CST_3; temp];
end
A_1 = [CST_3, zeros(size(CST_3,1),1)];
B_1 = [];
for i = 1:N
    B_1 = [B_1; G * [-wf; -mf]];
end


% w0 + w1 + w2 + ... <= fv_remain/delta
temp = eye(transNum);
for i = 1:N-1
    temp = [temp, eye(transNum)];
end
rows = size(temp,1);
A_22 = [zeros(rows, (N+1) * placeNum), temp, -1 * zeros(rows, 1)];
B_22 = fv_remain/delta - N * wf;



% w0 > zeta
%zeta = 0.001;
temp = zeros(1, (N+1) * placeNum + N * transNum);
temp(1, placeNum * (N+1) + 1 : placeNum * (N+1) + transNum) = ones(1, transNum);
A_3 = -1 * [temp, 1];
B_3 = -zeta;  % / delta;

% relaxation variable positive
A_4 = zeros(1, VarNum+1);
A_4(VarNum+1) = -1;
B_4 = 0;


% % sum(m_k+1) > epsilon
% temp = zeros(1, (N+1) * placeNum + N * transNum);
% temp(placeNum+1 : 2 * placeNum) = ones(1, placeNum);
% A_5 = -1* [temp, 1];
% B_5 = - epsilon;

% variable low bound
% temp = zeros(VarNum, VarNum+1);
% temp2 = zeros(VarNum,1);
% for i = placeNum*2+1 : placeNum*(N+1)
%     temp(i, i) = -1;
%     temp(i, VarNum+1) = -1;
% end
% for i = placeNum*(N+1)+transNum+1 : VarNum
%     temp(i, i) = -1;
%     temp(i, VarNum+1) = -1;
% end
% A_6 = temp;
% 
% start = placeNum*2 + 1;
% for i = 1 : N-1
%     for j = 1 : placeNum
%         temp2(start + (i-1)*placeNum +j) = -mf(j);
%     end
% end
% start = placeNum*(N+1) + transNum + 1;
% for i = 1:N-1
%     for j = 1 : transNum
%         temp2(start + (i-1)*transNum +j) = -wf(j);
%     end
% end
% B_6 = temp2;

lbVec  =  zeros(VarNum,1);
for i = 1 : N+1
    lbVec((i-1)*placeNum + 1 : i * placeNum) = -mf;
end
for i = 1 : N
    lbVec(placeNum * (N+1) + (i-1)*transNum + 1 : placeNum * (N+1) + i * transNum) = -wf;
end

temp = zeros(VarNum, VarNum+1);
for i = 1 : VarNum
    temp(i, i) = -1;
    temp(i, VarNum+1) = -1;
end
A_6 = temp;
B_6 = -1 * lbVec;


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
        lb(startIdx + transNum * (i-1) + j, 1) = -wf(j);
    end
end


% Up bound
%ub = Inf * ones(VarNum+1,1);

% %%%%%%%%%%%%%%%%%%%%%%%%%
% startIdx = 0;
% for i = 1 : N +1
%     if i > 2
%         for j = 1 : placeNum        
%             lb(startIdx + placeNum * (i-1) + j, 1) = -mf(j);
%         end
%     else
%         lb(startIdx + placeNum * (i-1) + j, 1) = -inf;
%     end
% end
% 
% startIdx = placeNum * (N+1);
% for i = 1 : N
%     if i > 1
%         for j = 1 : transNum
%             lb(startIdx + transNum * (i-1) + j, 1) = -wf(j) ;
%         end
%     else
%         lb(startIdx + transNum * (i-1) + j, 1) = -inf; 
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%

lb = [];
ub = [];

H = [];
% objective function
for i = 1 : N
    H = blkdiag(H, Q);
end
H = blkdiag(H, Z);
for i = 1 : N
    if i == 1
        H = blkdiag(H, -1 * R);
    else
        H = blkdiag(H, zeros(transNum));
    end
end


H = blkdiag(H, zeros(1));

f = zeros(VarNum + 1, 1);
f(VarNum + 1) = 1e10; %should be much bigger than Q,Z

Aeq = [Aeq_1; Aeq_2];
Beq = [Beq_1; Beq_2];

A = [A_1; A_22; A_3; A_4; A_6];
B = [B_1; B_22; B_3; B_4; B_6];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%using quadprog%%%%%%%%%%%%%
X0 = ones((N+1)*placeNum + N * transNum + 1, 1);
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
    fprintf('*');
    %m0 = m0
    %mf = mf
    %SimHPN_control_minFv(Post-Pre, m0, mf)     
    %halt

end

