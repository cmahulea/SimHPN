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

% SimHPN_control_onoff_plus() proportional ON-OFF controller, it is assumed that m0 > 0

%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre, Post: Pre- and post- incidence matrix
% m0:        Initial marking
% lambda:    Transitions firing rates
% delta:     Sampling period
% mf:        Final marking
% fv:        The firing count vector, if fv == [], using a minimal one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Output parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% s:  number of steps for reaching the final marking
% m:  |P|*(s+1) matrix, marking trajectory of places, m(:, k) is the
%     state in k_th step, m(:, 1) = m0, m(:, s+1) = mf
% u:  |T|*(s) matrix, control actions
% w:  |T|*(s) matrix, the controlled flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the transitions in a conflict relation are fired proportionally according
% to the firing count vector.

function [ s, m, u, w ] = SimHPN_control_onoff_plus( Pre, Post, m0, lambda, delta, mf, fv )

% Norm2 distance threshold to final state
qdis = 1e-6;

tStart = tic;

[placeNum, transNum] = size(Pre);
verified = 1;

% incidence matrix
C = Post - Pre;

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

%G = [DT -TA];

% get the starting time tick
tStart = tic;

if (size(fv, 1) == 0)
    % Compute minimal firing count vector
    fv = SimHPN_control_minFv(C, m0, mf);    
    fprintf('\nCorresponding minimal firing count vector is: %s\n ', mat2str(fv,4));  
else
    fprintf('\nUsing a given firing count vector: %s\n', mat2str(fv,4));
end

if isequal(fv, zeros(size(Pre,2), 1))
    s = 0;
    return;
end

% get the starting time tick
tStart = tic;

%%%%%%%%%%%%%%%%%%%Variables in LPP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = [w_{k}, m_{k+1}]'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = [ones(1, transNum), zeros(1, placeNum)];
s = 0;
m = m0;
u = [];
w = []; 
mpre = m0;
remainfv = fv;

% propositional firing for transitions in confilic
cm = SimHPN_control_cflBalance(Pre);

ctmax = SimHPN_control_conMatx(cm, fv);

% constraints m_k+1 - C * w_k * theta = m_k
A1 = [- C * delta, eye(placeNum)];

% w_k <= remain fv
A2 = [eye(transNum), zeros(transNum, placeNum)];

% flow constrains : DT * wk <= TA * mk
A3 = [DT, zeros(size(DT,1), placeNum)];    

% conflicting transitions are proportionally fired
if size(ctmax,1) ~= 0 
    A4 = [ctmax, zeros(size(ctmax,1), placeNum)];
else
    A4 = [];   
end

a = [A1;A2;A3;A4];

% bounds
lb = zeros(placeNum + transNum, 1);
ub = [];

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
    

while ( (mpre - mf)' * (mpre - mf) > qdis )
    s = s + 1;
    mk = mpre;

    % constrans value
    b1 = mk;
    b2 = remainfv / delta;
    b3 = TA*mk;
    b4 = zeros(size(A4,1), 1);

    b = [b1;b2;b3;b4];
        
    % invoke glpk
    [xopt, fmin, status, extra] = glpk (c, a, b, lb, ub, ctype, vartype, sense, param);

    mcur = xopt(transNum+1: transNum+placeNum, 1);
    wk = xopt(1:transNum, 1);

    m = [m, mcur];
    w = [w, wk];   
    u = [];
    
    mpre = mcur;
    remainfv = remainfv - wk * delta;
end

% get ending time tick
global tEscape;
tEscape = toc(tStart);
fprintf('Time consumed: %f millisecond, step: %d\n', tEscape * 1000, size(m,2) - 1);


