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

function SimHPN_Branch_Bound(ineq, pre, post, m0, lambda, AA, bb, ctype, vartype, eqs, control, contrtr, contrw, contrz)
%ineq = 1 --- minimizar
%ineq = 2 --- maximizar
global bound sol;

nps=size(pre,1);  % Number of places
nts=size(pre,2);  % Number of transitions

f = zeros(1,2*nts+nps);
if (ineq == 1)
    f(1) = 1;
else
    f(1) = -1;
end
if (control == 1)
    f(1:nts) = -contrw;
    f(nts+1:nts+nps) = contrz;
end

% add the equality corresponding to the place that limits the flow of the
% transitions
for i = 1 : length(eqs)
    add = eqs{i};
    temp = zeros(1,2*nts+nps);
    temp(add(1)) = 1;
    temp(nts+add(2)) = -lambda(add(1))/pre(add(2),add(1));
    AA = [AA;temp];
    bb = [bb;0];
    ctype = sprintf('%sS',ctype);
end

  
[X,dummy,status] = glpk(f,AA, bb, zeros(size(AA,2),1), [], ctype, vartype, 1);

if ((status ~= 5) || (X(1) >= bound)) && (ineq == 1)
    return;
end
if ((status ~= 5) || (X(1) <= bound)) && (ineq == 2)
    return;
end

% if ((EXITFLAG <= 0) || (X(1) >= bound)) && (ineq == 1)
%     return;
% end
% if ((EXITFLAG <= 0) || (X(1) <= bound)) && (ineq == 2)
%     return;
% end

fx = X(1:nts);
mx = X(nts+1:nts+nps);

%compute set nt - not satisfied transitions
Pre1 = zeros(size(pre,1),size(pre,2));
for i = 1 : size(pre,2)
    inputP = find(pre(:,i));
    temp = mx(inputP);
    for j = 1 : length(inputP)
        temp(j) = temp(j)/pre(inputP(j),i);
    end
    [mi,ii] = min(temp);
    Pre1(inputP(ii),i) = 1/pre(inputP(ii),i);
end
Pre1 = Pre1';
fxx = diag(lambda)*Pre1*mx;
nt = find(abs(fxx-fx) > 1e-3);

%remove from nt all controllable transitions
if ((ineq == 2) && (control == 1))
    nt = setdiff(nt,contrtr);
end

if isempty(nt)
    if (ineq == 1)
        bound = min(bound,X(1));
        sol = X;
    else
        bound = max(bound,X(1));
        sol = X;
    end
    return;
end

%take a t from nt
t = nt(1);
inputPlaces = find(pre(:,t) > 0);
eqs{length(eqs)+1} = 1;
for i = 1 : length(inputPlaces)
    eqs{length(eqs)} = [t inputPlaces(i)];
    SimHPN_Branch_Bound(ineq, pre, post, m0, lambda, AA, bb, ctype, vartype, eqs, control, contrtr, contrw, contrz);
end
