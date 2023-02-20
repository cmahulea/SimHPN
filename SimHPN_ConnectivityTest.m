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


function [G,Keq,flag] = SimHPN_ConnectivityTest(Pre,Post,Tc,L,m0)

%Reducing the configurations to analyze
%1) Removing implicit arcs
imp_arc = SimHPN_ImplicitArcs(Pre,Post,m0);%;

Pre_red = Pre;
if ~isempty(imp_arc)
    for h=1:size(imp_arc,1)
        Pre_red(imp_arc(h,1),imp_arc(h,2))=0;
    end
end

conf_red = Configurations(Pre_red);%;   %Set of reduced configurations.
                                     %Only the ones that don't contain
                                     %implicit arcs are preserved    


                                     
%2)Preserving only the configurations in which E*~=0   

[Keq, A, Aeq] = SimHPN_ConfigurationsWEquilibria(Pre,Post,m0,L,Tc,conf_red);

[V,E] = SimHPN_Connectivity(Pre,Post,m0,L,Tc,Keq);


%3) Verifying Connectivity over the reduced set of configurations

s = E(:,1)';
t = E(:,2)';
G = graph(s,t);
flag = 0;
bins = conncomp(G);
if size(unique(bins),2)==1
    flag=1;
end
    
return

