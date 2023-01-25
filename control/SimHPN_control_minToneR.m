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

% compute the minimal time for reaching mf from m0, in the case that
% mf and m0 are in the same region and the flow of transitions are CONSTANT

%%%%%%%%%%%%%%%%%%%%%%%%%Output Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  minT:   The minimal time reaching mf from m0
%  w   :   the controlled firing flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [minT, w] = SimHPN_control_minToneR(Pre, Post, lambda, m0, mf)

[placeNum, transNum] = size(Pre);


cfgM0 = SimHPN_control_CfgMatrixGen(Pre, Post, m0);
cfgMf = SimHPN_control_CfgMatrixGen(Pre, Post, mf);


A = [Post - Pre, zeros(placeNum, 1);
     eye(transNum), -1 * lambda .* (cfgM0 * m0);
     eye(transNum), -1 * lambda .* (cfgM0 * mf)];


 % the sense of each constraint in the constraint matrix
for i = 1 : (2 * transNum + placeNum)
   if i <= placeNum
       ctype(i, 1) = 'S';
   else
       ctype(i,1) = 'U';
   end
end

% types of the variables
for i = 1 : transNum + 1
    vartype(i,1) = 'C';
end

% objective function coefficients
co = zeros(1, transNum + 1);
co(1, transNum + 1) = 1;

lb = zeros(transNum + 1, 1);
ub = [];

for i = 1 : (placeNum + 2 * transNum)
    if i <= placeNum
        b(i, 1) = mf(i) - m0(i);
    else
        b(i, 1) = 0;
    end
end


param.msglev = 0;
[xopt, fopt, status, extra] = glpk (co, A, b, lb, ub, ctype, vartype, 1, param);

if status == 5
    minT = fopt;
    w(:, 1) = xopt(1:transNum);
else
    minT = 0
    w(:, 1) = zeros(transNum, 1)
    fprintf('LPP doesnot find solution!\n');
    halt
end






