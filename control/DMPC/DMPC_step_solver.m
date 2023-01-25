% solve the optimization problem %    This is part of SimHPN Toolbox, for Matlab 2010b or newer.
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

of subsystems in each step
% Input parameter: subsys, lambda, N, Q, Z, phase
%                  phase = 1 normal procedure
%                  phase = 2 when some buffer places are emptied
%                  alpha, 0 <= alpha < 1
%                  eps: a small positive number
% Output parameter: subsys_updated (the states of buffers are not updated)
%                   subsys_updated{i}.Pre the Pre matrix of subsystem i 
%                   subsys_updated{i}.Post the Pre matrix of subsystem i 
%                   subsys_updated{i}.p the places in subsytem i
%                   subsys_updated{i}.t the transitios in subsytem i
%                   subsys_updated{i}.it the interface transitios in subsytem i
%                   subsys_updated{i}.b the buffer places in subsytem i
%                   subsys_updated{i}.bm the state of buffer places
%                   subsys_updated{i}.bmf the final state of buffer places 
%                   subsys_updated{i}.m the state of subsystem i 
%                   subsys_updated{i}.mf the final state of subsystem i 
%                   subsys_updated{i}.lambda the firing rate of subsystem i
%                   subsys_updated{i}.w the flow of subsystem i


function subsys_updated = DMPC_step_solver(subsys, delta, N, Q, Z, phase, alpha, eps)

global Pre_G;
global Post_G;

subsys_updated = subsys;

[placeNum, transNum] = size(subsys.Pre);
m0 = subsys.m;
mf = subsys.mf;
wf = zeros(transNum,1);

% G: the flow constraint
DT = [];
TA = [];
for (i = 1:placeNum)
    for (j = 1:transNum)
        if (subsys.Pre(i, j) > 0)
            tmpD = zeros(1, transNum);
            tmpD(j) = 1;
            tmpTA = zeros(1, placeNum);
            tmpTA(i) = subsys.lambda(j)/subsys.Pre(i, j);
            DT = [DT;tmpD];
            TA = [TA;tmpTA];
        end
    end
end
G = [DT -TA];

C = subsys.Post - subsys.Pre;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X = [m(k)-mf, m(k+1)-mf,..., m(k+N)-mf, w(k)-wf, w(k+1)-wf,..., w(k+N-1)-wf]'
% k = 0, 1, 2...
% wf = 0
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

Aeq_1 = [CST_1, -1 * CST_2];
Beq_1 = [];
for i = 1 : N
    Beq_1 = [Beq_1; C * delta * wf];
end

% constrains matrix for: m(k)-mf = m0 - mf
Aeq_2 = [eye(placeNum), zeros(placeNum,N*placeNum + N*transNum)];
Beq_2 = m0 - mf;


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

relax = 1e-10;

% mf(pi) >= mo(pi) --> m_k+1(pi) <= mf(pi); otherwise --> m_k+1(pi) >= mf(pi)
temp = zeros(placeNum, placeNum);
for i = 1 : placeNum
    if m0(i,1) <= mf(i,1)
        temp(i,i) = 1;
    else
        temp(i,i) = -1;
    end
end
A_2 = zeros(placeNum, placeNum);
for i = 1 : N
    A_2 = blkdiag(A_2 , temp);
end
A_2 = [A_2, zeros(size(A_2,1), N * transNum)];
B_2 = zeros(size(A_2,1), 1);

% mf(pi) >= mo(pi) --> m_k+1(pi) >= m_k(pi); otherwise --> m_k+1(pi) <= m_k(pi)
temp = zeros(placeNum, placeNum);
for i = 1 : placeNum
    if m0(i,1) <= mf(i,1)
        temp(i,i) = 1;
    else
        temp(i,i) = -1;
    end
end
temp2 = -temp;

Part_1 = [];
for i = 1 : N
    Part_1 = blkdiag(Part_1, temp);
end
Part_2 = [];
for i = 1 : N
    Part_2 = blkdiag(Part_2, temp2);
end

A_3 = [[Part_1, zeros(size(Part_1,1), placeNum)] + [zeros(size(Part_1,1), placeNum), Part_2], ... 
       zeros(size(Part_1,1), N * transNum)];
B_3 = zeros(size(A_3,1), 1) +  relax * ones(size(A_3,1), 1);


% consider the constraint from the buffer place in the first step
A_4 = [];
B_4 = [];

for biSub = 1 : size(subsys.b,2)
    biGlb = subsys.b(1, biSub);  
    
    for j = 1 : size(subsys.it,2)
        itGlb = subsys.it(j);
        
        if Pre_G(biGlb, itGlb) > 0
            newRow = zeros(1, transNum);
            itSub = find(subsys.t == itGlb);
            newRow(1,itSub) = 1;
            A_4 = [A_4; newRow];
            fmax = subsys.lambda(itSub) * subsys.bm(biSub) / Pre_G(biGlb, itGlb);
            B_4 = [B_4; fmax]; 
        end
    end
end
rowN = size(A_4,1);
A_4 = [zeros(rowN, placeNum * (N+1)) A_4 zeros(rowN, transNum * (N-1))];


% in phase 2, do not empty any emptied buffer place
A_6 = [];
B_6 = [];

if phase == 2
    for biSub = 1 : size(subsys.b,2)
        biGlb = subsys.b(1, biSub);  
        
        for j = 1 : size(subsys.it,2)
            itGlb = subsys.it(j);
            newRow = zeros(1, transNum);
            if Pre_G(biGlb, itGlb) > 0
                itSub = find(subsys.t == itGlb);
                newRow(1,itSub) = delta * Pre_G(biGlb, itGlb);
                %itGlb = itGlb
                %biGlb = biGlb
            end
        end
        if ~isequal(newRow, zeros(1, transNum))
            A_6 = [A_6; newRow];
            B_6 = [B_6; alpha * subsys.bm(biSub)];
        end
    end
    rowN = size(A_6,1);
    A_6 = [zeros(rowN, placeNum * (N+1)) A_6 zeros(rowN, transNum * (N-1))];
    
end


% Low bound 
startIdx = 0;
for i = 1 : N +1
    for j = 1 : placeNum        
        lb(startIdx + placeNum * (i-1) + j, 1) = -mf(j) - relax;
    end
end

startIdx = placeNum * (N+1);
for i = 1 : N
    for j = 1 : transNum
        lb(startIdx + transNum * (i-1) + j, 1) = -wf(j)  - relax;
    end
end

% upper bound
ub = Inf * ones(VarNum,1);

H = [];
% objective function
for i = 1 : N
    H = blkdiag(H, Q);
end
H = blkdiag(H, Z);

R = zeros(transNum);
for i = 1 : N
    H = blkdiag(H, R);
end

A_5 = [];
B_5 = [];
f = zeros(VarNum, 1);
temp = [];

% consdier the case when some of its output buffers are (almost) emptied
if phase == 2  
    %fprintf('--------------------entered phase 2----------------------\n');
    startIdx = placeNum * (N+1);
    
    for biSub = 1 : size(subsys.b,2)        
        if subsys.bm(biSub) <= eps %subsys.bmf(biSub) % eps
            biGlb = subsys.b(1, biSub);
            
            for j = 1 : size(subsys.it,2)
                itGlb = subsys.it(j);
                if Post_G(biGlb, itGlb) > 0 % output buffer
                    % find the index of tj in subsys.t  
                    itSub = find(subsys.t == itGlb);
                    f(startIdx + itSub) = -1;                    
                end
            end
        end
    end   
    H = zeros(size(H));
else
    %fprintf('--------------------entered phase 1----------------------\n');
end


Aeq = [Aeq_1; Aeq_2];
Beq = [Beq_1; Beq_2];

A = [A_1; A_2; A_3; A_4; A_5; A_6];
B = [B_1; B_2; B_3; B_4; B_5; B_6];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%using quadprog%%%%%%%%%%%%%
%X0 = rand((N+1)*placeNum + N * transNum, 1);
%X0(1:placeNum,1) = m0 - mf;
X0 = [];
options = optimset('LargeScale', 'OFF', 'Display', 'OFF');
%options=optimset('TolCon', 1e-004);
% sH = size(H)
% sf = size(f)
% sA = size(A)
% sB = size(B)
% sAeq = size(Aeq)
% sBeq = size(Beq)
% sLb = size(lb)
% sUb = size(ub)

global FVAL;

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

%subsys_updated.m = [subsys_updated.m, mArray(:, 1)];
%subsys_updated.w = [subsys_updated.w, wArray(:, 1)];

subsys_updated.m = mArray(:, 1);
subsys_updated.w = wArray(:, 1);


%subsys_updated.m
%wwwww = subsys_updated.w

%mArray = mArray
%wArray = wArray

%mArray = mArray
%wArray = wArray  

if EXITFLAG == -2
    fprintf('***********************************************************\n');
    %m0 = m0
    %mf = mf
    %SimHPN_control_minFv(Post-Pre, m0, mf)     
    halt
else
    %fprintf('EXITFLAG = %s \n', num2str(EXITFLAG));
end
