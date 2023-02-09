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

function [Pi,Ti,Flag] = SimHPN_Influ(Pre,Post,Tc)

numT = size(Pre,2);
numP = size(Pre,1);
C = Post-Pre;
Flag = 0;
if isempty(Tc)
    Pi=[];Ti=[];Flag=0;
    disp('The set of controllable transitions is empty. Hence, no node is influenced');
elseif max(Tc)>numT
    disp('There is an error in Tc');
    Pi=[];Ti=[];Flag=0;
else
    Ti = [Tc];
    A=[];B=[];
    for k=1:size(Tc,2)
        A = [A find(C(:,Tc(k)))']  ;     %Find all places that are input or output to a controllable transition (no balanced self loops)
    end
    A = unique(A);
    sl_p = find(all(C==0,2))';          %Find all of the balanced self loops 
    if ~isempty(sl_p)
        for j=1:size(sl_p,2)
            pre_slp = find(Pre(sl_p(j),:))';
            if all(ismember(pre_slp,Tc))
                B = [B sl_p(j)];
            end
        end
    end
    Pi = [A B];
    rp_f = 0;
    while rp_f==0
        Pa = Pi;
        Ta = Ti;
        Tn=[];Pn=[];
        for i=1:size(Pa,2)
            ptp = setdiff(find(Pre(Pa(i),:)),Ta);
            if ~isempty(ptp)
                for j=1:size(ptp)
                    if all(ismember(find(Pre(:,ptp(j)))',Pa))
                        Tn = [Tn ptp(j)];
                    end
                end
            end
        end
        Tn = unique(Tn);
        if ~isempty(Tn)
            Ti = unique([Ta Tn]);
            for k=1:size(Tn,2)
                pretn = find(Pre(:,Tn(k)))';
                postn = find(Post(:,Tn(k)))';
                prepostn = unique([pretn postn]);
                Pn = unique([Pn prepostn]);
            end
            Pi = unique([Pa Pn]);
        end
        if all(ismember(Pi,Pa))
            rp_f = 1;
        end  
    end
end

if all(ismember([1:numP],Pi))
    Flag = 1;
end

disp('The set of influenced nodes are')
Pi
Ti
if Flag==1
    disp('Influence is total')
else
    disp('Influence is not total, thus, the net is not NRC')
end

return

