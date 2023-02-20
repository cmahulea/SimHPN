%Given a TCPN system (N,L,m0), this algorithm computes the set of all the
%configurations of the system s.t. the corresponding Ei* is not an emptyset.
function [Keq, A, Aeq] = SimHPN_ConfigurationsWEquilibria(Pre,Post,mu0,L,Tc,conf)

np = size(Pre,1);
nt = size(Pre,2);
C = Post-Pre;
nBy = np-rank(C);
nTc = size(Tc,2);
Tnc = setdiff(1:nt,Tc); nTnc = size(Tnc,2);
Tr = eye(nt);
TcM = Tr(Tc,:);
TncM = Tr(Tnc,:);
%conf = Configurations(Pre);

By = null(C','r');

num_conf = size(conf,1);
Req = [];
T_sync = [];
for k=1:nt
    if nnz(Pre(:,k))>1
        T_sync = [T_sync k];
    end
end
%%
for j=1:num_conf
%Que configuración quieres analizar?
conf_i = conf(j,:);
Pi_i = compPi(Pre,conf_i);

%Matriz for the equalities,... Class(m0), equilibrium and 
%control of uncontrollables equal to 0
Aeq =[By' zeros(nBy,nt) zeros(nBy,1); 
             C*L*Pi_i -C zeros(np,1); 
             zeros(nTnc,np) TncM zeros(nTnc,1)];
beq = [By'*mu0;
       zeros(np,1);
       zeros(nTnc,1)];
%Ri corresponds to the matrix defining the inequalities of the
%region
Ri=[];
for t=1:size(T_sync,2)
    P_in = find(Pre(:,T_sync(t)))';
    for s=1:size(P_in,2)
        if conf_i(T_sync(t))~=P_in(s)
            v = zeros(1,np);
            v(conf_i(T_sync(t)))=1/Pre(conf_i(T_sync(t)),T_sync(t));
            v(P_in(s))=-1/Pre(P_in(s),T_sync(t));
            Ri = [Ri;v];
        end
    end
end
r = size(Ri,1);

A = [Ri zeros(r,nt) zeros(r,1);
    -L*Pi_i eye(nt) ones(nt,1);
    zeros(nTc,np) -TcM ones(nTc,1)];
b = zeros(size(A,1),1);

f = [zeros(1,np+nt) -1];
options = optimoptions('linprog','Display','none');
x = linprog(f,A,b,Aeq,beq,zeros(np+nt+1,1),150*ones(np+nt+1,1),options);

if isempty(x)
    x = ones(np+nt+1,1);
end
if x(np+nt+1)>0
    if (abs(beq - Aeq*x)<=.001*ones(size(Aeq,1),1))
        Req = [Req;j];
    end
end
end

Keq = conf(Req,:);