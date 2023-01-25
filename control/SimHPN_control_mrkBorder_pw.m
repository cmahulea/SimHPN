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

% find the marking on the board, the flow of transitions are constant
% within a region

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mb_pw: placeNum * S matrix, S is the number of borders being crossed
%        each colum is a marking on the board in order
% minTmb_pw : 1 * (S + 1) matirx 
%             the minimal time for reaching the marking on boarder from
%             the last marking. the first marking is m0, the last marking
%             is mf
% w: transNum * (S + 1) matrix 
%    controlled flow (CONSTANT) between states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X = [m0 m1 ... mf(= m_s+1) x1 x2 ... x_s+1 t1 t2 ... t_s+1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mb_pw, minTmb_pw, w] =  SimHPN_control_mrkBorder_pw(Pre, Post, lambda, m0, mf)

[placeNum, transNum] = size(Pre);

% find the border being corossed to reach mf from m0
% mb is the makring on border in the case of constant flow

[mb, PiRespRowDif] = SimHPN_control_mrkBorder(Pre, Post, m0, mf);
S = size(mb, 2); % number of border being crossed

%mb = mb
%NumOfBorderCrossed = S

if (S == 0)
    mb_pw = [];
    minTmb_pw = [];
    w = [];
    return;
end

C = Post - Pre;
CC = [];
for (si = 1 : (S+1))
    CC = blkdiag(CC, C);
end

allOneDiag = [];
allMinOneDiag = [];

for (si = 1 : (S+1))
    allOneDiag = blkdiag(allOneDiag, eye(placeNum));
    allMinOneDiag = blkdiag(allMinOneDiag, -1 * eye(placeNum));
end

BB = [allOneDiag,  zeros(size(allOneDiag,1), placeNum)] + [zeros(size(allMinOneDiag, 1), placeNum), allMinOneDiag ];

DD = [];
for (si = 1 : S)
    
    if (si == S)
        blk = SimHPN_control_CfgMatrixGen(Pre, Post, mf) - SimHPN_control_CfgMatrixGen(Pre, Post, mb(:, S));   
    else
        blk = SimHPN_control_CfgMatrixGen(Pre, Post, mb(:, si + 1)) - SimHPN_control_CfgMatrixGen(Pre, Post, mb(:, si));
    end  
    
    DD = blkdiag(DD, blk);
end

DDr = size(DD, 1);
DD = [zeros(DDr, placeNum), DD, zeros(DDr, placeNum)];

for (i = 1 : placeNum)
    if (m0(i) < mf(i))
        blkA(i, i) = 1;
    else
        blkA(i, i) = -1;
    end
end

blkB = -1 * blkA;

tmpA = []; tmpB = [];

for(si = 1 : (S+1))
    tmpA = blkdiag(tmpA, blkA);
    tmpB = blkdiag(tmpB, blkB);
end

EE = [tmpA, zeros(size(tmpA, 1), placeNum)] + [zeros(size(tmpB), placeNum), tmpB];

% m0 = m_0, mf = m_(s+1)    
FF = [eye(placeNum), zeros(placeNum, (S+1) * placeNum) ; zeros(placeNum, placeNum * (S + 1)), eye(placeNum)];

Aeq = [BB, CC, zeros((S+1) * placeNum, S+1); FF, zeros(2 * placeNum, (transNum + 1) * (S + 1)) ; DD, zeros(S * transNum, (transNum + 1) * (S + 1))];
Beq = [zeros((S+1) * placeNum ,1); m0; mf; zeros(S * transNum, 1)];

A = [EE, zeros((S+1)*placeNum, (transNum + 1) * (S +1))];
B = zeros((S+1)*placeNum, 1);


X0 = zeros(placeNum * (S+2) + transNum * (S+1) + S + 1, 1);
for (si= 1:S+2)
    X0((si-1) * placeNum + 1 : placeNum * (si)) = m0;
end


LB = zeros(placeNum * (S+2) + transNum * (S+1) + S + 1, 1);
UB = [];

options=optimset('Algorithm','active-set', 'Display', 'off');
[X, fval, exitflag] = fmincon(@(X) objFunc(X, Pre, S),X0,A,B,Aeq,Beq,LB,UB, @(X) nonLinFunc2(X, Pre, Post, S, m0, mf, mb, lambda), options);

if (exitflag == -2)
    fprintf('PP no solution\n');
    exitflag = exitflag
    halt;
end

mb_pw = zeros(placeNum, S);
for (si = 1:S)
	mb_pw(1 : placeNum, si) = X( si * placeNum + 1 : (si+1) * placeNum);
end

minTmb_pw = X((S+2) * placeNum + (S+1) * transNum + 1 : (S+2) * placeNum + (S+1) * transNum + S + 1)';

for (si = 1 : S+1)
    %w((si-1)*transNum + 1 : si*transNum, si) = X((S+2) * placeNum + (si-1)*transNum + 1 : (S+2) * placeNum + si*transNum);
    w(1 : transNum, si) = X((S+2) * placeNum + (si-1)*transNum + 1 : (S+2) * placeNum + si*transNum);
end

% objection function
function f = objFunc(X, Pre, S)
[placeNum, transNum] = size(Pre);
f = [zeros(1, (S+2)*placeNum + (S+1) * transNum), ones(1, (S+1))] * X;


%non-linear constrain function
function [c,ceq] = nonLinFunc(X, Pre, Post, S, m0, mf, mb, lambda)

ceq = [];
[placeNum, transNum] = size(Pre);

blk = eye(transNum);
st1 = [];
for (si = 1: (S+1))
    st1 = blkdiag(st1, blk);
end

st2 = [];
temp = [];

for (si = 1 : (S+1))
    st2((si - 1) * transNum + 1 : (si * transNum), 1) = lambda;
end

st3 =[];
temp = [];

for (si = 0 : (S))
    if (si == 0)
        cfg = SimHPN_control_CfgMatrixGen(Pre, Post, m0);
    %elseif (si == (S + 1) )
    %    cfg = SimHPN_control_CfgMatrixGen(Pre, Post, mf);
    else
        cfg = SimHPN_control_CfgMatrixGen(Pre, Post, mb(:, si));
    end
    temp = blkdiag(temp, cfg);
end

st3 = [temp, zeros(size(temp, 1), placeNum)];
st4 = [zeros(size(temp, 1), placeNum), temp];


MV = X(1:(S+2)*placeNum, 1);
XV = X((S+2)*placeNum + 1 : ((S+2)*placeNum + (S+1)*transNum), 1);
TV = X(((S+2)*placeNum + (S+1)*transNum) + 1 :((S+2)*placeNum + (S+1)*(transNum + 1)), 1); 

exp1 = st2 .* (st3 * MV);
exp2 = st2 .* (st4 * MV);


for (si = 1:S+1)
    idx = (si - 1) * transNum;
    for (j = 1:transNum)
        temp1 = exp1(idx + j, 1) * TV(si);
        temp2 = exp2(idx + j, 1) * TV(si);
        
        exp1(idx + j, 1) = temp1;
        exp2(idx + j, 1) = temp2;        
    end
end


c1 = st1 * XV - exp1;
c2 = st1 * XV - exp2;


c = [c1 ; c2];


% X = [m0 m1 ... mf(= m_s+1) x1 x2 ... x_s+1 t1 t2 ... t_s+1]
%non-linear constrain function
function [c,ceq] = nonLinFunc2(X, Pre, Post, S, m0, mf, mb, lambda)

ceq = [];
c = [];
[placeNum, transNum] = size(Pre);

for (k = 1: S+1)
    startIdx = (S+2)*placeNum + (k - 1) * transNum;
    XVk = X(startIdx +1 : startIdx + transNum, 1);
    
    startIdx = (S+2)*placeNum + (S+1)*transNum;
    TVk = X(startIdx + k, 1);
    
    startIdx = (k-1)*placeNum;
    MVk_1 = X(startIdx + 1 : startIdx + placeNum, 1);
    MVk = X(startIdx + placeNum + 1 : startIdx + 2 * placeNum, 1);
  
    cfgk_1 = SimHPN_control_CfgMatrixGen(Pre, Post, MVk_1);
    cfgk = SimHPN_control_CfgMatrixGen(Pre, Post, MVk);
    c = [XVk - lambda .* (cfgk_1 * MVk_1) * TVk;
         XVk - lambda .* (cfgk * MVk) * TVk;
         c];
end




