function [V,E] = SimHPN_Connectivity(Pre,Post,mu0,L,Tc,Keq)

sK = size(Keq,1);
V=1:sK;   %Vertices del grafo
E = [];
C = Post-Pre;
By = null(C','r');
for i=1:sK-1
    for j=i+1:sK
        [Ctv] = SimHPN_EiEjConnected(Pre,Post,By,mu0,L,Tc,Keq(i,:),Keq(j,:));
        if Ctv == 1
            E = [E;[i,j]];
        end
    end
end
