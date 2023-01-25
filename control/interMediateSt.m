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

function m = interMediateSt( Pre, Post, m0, mf, lambda, delta, fv )

[placeNum, transNum] = size(Pre);
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

G = [DT -TA];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for model 3.
% x = [m, sigma, w]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m = m0 + C * sigma
A1 = [eye(placeNum), -C, zeros(placeNum,transNum)];
b1 = m0;

% sigma <= fv
A2 = [zeros(transNum, placeNum), eye(transNum), zeros(transNum, transNum)];
b2 = fv;

% flow constrains :  TA * m - DT * w >= 0
A3 = [TA, zeros(size(TA,1), transNum), -DT]; 
b3 = zeros(size(A3,1),1);

% proportionally fired
cm = SimHPN_control_cflBalance(Pre);
ctmax = SimHPN_control_conMatx(cm, fv);

A4 = [zeros(size(ctmax,1), placeNum), zeros(size(ctmax,1), transNum), ctmax];
b4 = zeros(size(A4,1),1);

% constraint matrix
A = [A1;A2;A3;A4];
b = [b1;b2;b3;b4];

% bounds
lb = zeros(placeNum + 2 * transNum,1);
%lb(1:placeNum,1) = 1e-5 * ones(placeNum,1); % not negative

ub = inf * ones(placeNum + 2 * transNum,1);
ub(placeNum + 3) = 0;
ub(placeNum + 4) = 0;

% type of variables
for i = 1 : placeNum + 2 * transNum
    vartype(i,1) = 'C';
end

for i = 1 : size(A1,1)
    ctype1(i,1) = 'S';
end

for i = 1 : size(A2,1)
    ctype2(i,1) = 'U';
end

for i = 1 : size(A3,1)
    ctype3(i,1) = 'L';
end

for i = 1 : size(A4,1)
    ctype4(i,1) = 'S';
end

ctype = [ctype1; ctype2; ctype3; ctype4];

c = zeros(1, placeNum + 2 * transNum);

%c(placeNum+ transNum + 1) = 1;
%c(placeNum+ transNum + 2) = 1;
%c(placeNum+ transNum + 3) = 1;
%c(placeNum+ transNum + 4) = 1;

c(placeNum+ transNum + 1: placeNum+ 2* transNum) = 1;

% sense, maximize
sense = -1; 
param.msglev = 1;


[xopt, fmin, status, extra] = glpk (c, A, b, lb, ub, ctype, vartype, sense, param);

%xopt

mT = xopt(1:placeNum)
sigT = xopt(placeNum+1:placeNum+transNum)
wT = xopt(placeNum+transNum + 1:placeNum + 2 * transNum)


SimHPN_control_onoff_plus(Pre, Post, m0, lambda, delta, mT, sigT );

SimHPN_control_onoff_plus(Pre, Post, mT, lambda, delta, mf, fv - sigT );


%SimHPN_control_onoff_plus(Pre, Post, m0, lambda, delta, mf, fv );




    
