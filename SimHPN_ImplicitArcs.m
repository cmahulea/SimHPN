%    This is part of SimHPN Toolbox, for Matlab 2010b or newer.
%
%    Copyright (C) 2023 SimHPN developing team. For people, details and citing
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
function imp_arc = SimHPN_ImplicitArcs(Pre,Post,mu0)

np = size(Pre,1);
nt = size(Pre,2);
C = Post-Pre;
P = [1:np];

T_sync = [];
for k=1:nt
    if nnz(Pre(:,k))>1
        T_sync = [T_sync k];
    end
end

%que arco quieres analizar?? uno a uno con las transiciones de T_sync
imp_arc = [];
for t=1:size(T_sync,2)
    ct = T_sync(t);
    P_ = find(Pre(:,ct))';
    for p=1:size(P_,2)
        cp = P_(p);
        Pp = setdiff(P,cp);
        A = [C(Pp,:)';-Pre(Pp,t)'];
        b = [C(cp,:)';-Pre(cp,ct)'];
        f = mu0(Pp);
        Aeq = [];
        beq = [];
        options = optimoptions('linprog','Display','none');
        [y,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,zeros(np-1,1),[],options);
        if exitflag==1
            if mu0(cp)>=f'*y
                imp_arc = [imp_arc;cp ct];
            end
        end
    end
end

return;