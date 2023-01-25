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

% compute the intermediate state between m0 and mf to minimize the time reaching mf
% m0 and mf are in the same region, and the flows bewteen states are constant here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% md: the intermediate state added between m0 and mf
% t1: the minimal time reaching md from m0
% t2: the minimal time reaching mf from md
% x1: the accumulated controlled flow reaching md from m0
% x2: the accumulated controlled flow reaching mf from md
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% X = [m0, md, mf, x1, x2, t1, t2]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [md, t1, t2, x1, x2] = SimHPN_control_interState_pw(Pre, Post, lambda, m0, mf)

[placeNum, transNum] = size(Pre);

% find the border being corossed to reach mf from m0

S = 1; % only add one intermediate state

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

BB = [allOneDiag,  zeros(placeNum * (S + 1), placeNum)] + [zeros(placeNum * (S + 1), placeNum), allMinOneDiag ];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

EE = [tmpA, zeros(placeNum * (S + 1), placeNum)] + [zeros(placeNum * (S + 1), placeNum), tmpB];


% m0 = m_0, mf = m_(s+1)    
FF = [eye(placeNum), zeros(placeNum, (S+1) * placeNum) ; zeros(placeNum, placeNum * (S + 1)), eye(placeNum)];


Aeq = [BB, CC, zeros((S+1) * placeNum, S+1); FF, zeros(2 * placeNum, (transNum + 1) * (S + 1))];
Beq = [zeros((S+1) * placeNum ,1); m0; mf];


A = [EE, zeros((S+1)*placeNum, (transNum + 1) * (S +1))];
B = zeros((S+1)*placeNum, 1);


X0 = ones(placeNum * (S+2) + transNum * (S+1) + S + 1, 1);

X0(1 : placeNum, 1) = m0;
X0(placeNum + 1 : 2 * placeNum,1) = mf; % + rand * (mf - m0);
X0(2 * placeNum + 1 : 3 * placeNum, 1) = mf;



LB = zeros(placeNum * (S+2) + transNum * (S+1) + S + 1, 1);
UB = [];

options=optimset('Algorithm','active-set', 'TolCon', 1e-005, 'Display', 'off');
[X, fval, exitflag] = fmincon(@(X) objFunc(X, Pre, S),X0,A,B,Aeq,Beq,LB,UB, @(X) nonLinFunc(X, Pre, Post, S, m0, mf, lambda), options);

if (exitflag == -2)
    fprintf('PP no solution\n');
    m0 = m0
    mf = mf
    exitflag = exitflag
    %halt;
end

md = X(placeNum + 1 : 2 * placeNum);
x1 = X(placeNum * 3 + 1 : placeNum * 3 + transNum);
x2 = X(placeNum * 3 + transNum + 1 : placeNum * 3 + 2 * transNum);
t1 = X(3*placeNum + 2 * transNum + 1);
t2 = X(3*placeNum + 2 * transNum + 2);


% objection function
function f = objFunc(X, Pre, S)
[placeNum, transNum] = size(Pre);
f = [zeros(1, (S+2)*placeNum + (S+1) * transNum), ones(1, (S+1))] * X;


%non-linear constrain function
function [c,ceq] = nonLinFunc(X, Pre, Post, S, m0, mf, lambda)

ceq = [];
[placeNum, transNum] = size(Pre);

MV = X(1:(S+2)*placeNum, 1);
XV = X((S+2)*placeNum + 1 : ((S+2)*placeNum + (S+1)*transNum), 1);
TV = X(((S+2)*placeNum + (S+1)*transNum) + 1 :((S+2)*placeNum + (S+1)*(transNum + 1)), 1); 

MV0 = MV(1:placeNum, 1);
MVd = MV(placeNum + 1 : 2*placeNum, 1);
MVf = MV(2*placeNum +1 : 3*placeNum, 1);

cfg = SimHPN_control_CfgMatrixGen(Pre, Post, m0);

XV1 = XV(1:transNum, 1);
XV2 = XV(transNum + 1 : 2 * transNum, 1);

TV1 = TV(1, 1);
TV2 = TV(2, 1);

c1 = XV1 - lambda .* (cfg * m0) * TV1;
c2 = XV1 - lambda .* (cfg * MVd) * TV1;
c3 = XV2 - lambda .* (cfg * MVd) * TV2;
c4 = XV2 - lambda .* (cfg * mf) * TV2;

%MstartIdx = 0;
%XstartIdx = (S+2) * placeNum;
%TstartIdx = (S+2) * placeNum + (S+1) * transNum;

%c1 = X(XstartIdx + 1: XstartIdx + transNum, 1) - lambda .* (cfg * m0) * X(TstartIdx + 1, 1);
%c2 = X(XstartIdx + 1: XstartIdx + transNum, 1) - lambda .* (cfg * X(MstartIdx + placeNum + 1 : MstartIdx + 2 *placeNum)) * X(TstartIdx + 1, 1);
%c3 = X(XstartIdx + transNum + 1: XstartIdx + 2 * transNum, 1) - lambda .* (cfg * X(MstartIdx + placeNum + 1 : MstartIdx + 2 *placeNum)) * X(TstartIdx + 2, 1);
%c4 = X(XstartIdx + transNum + 1: XstartIdx + 2 * transNum, 1) - lambda .* (cfg * mf) * X(TstartIdx + 2, 1);

c = [c1 ; c2; c3; c4];


