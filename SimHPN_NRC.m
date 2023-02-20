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

function [Flag,sc1,sc2,sc3,Pi,Ti] = SimHPN_NRC(Pre,Post,Tc)

numT = size(Pre,2);
numP = size(Pre,1);
C = Post-Pre;
Flag = 0;
if isempty(Tc)
    Pi=[];Ti=[];Flag=0;
    disp('The set of controllable transitions is empty. Hence, the net is not NRC');
elseif max(Tc)>numT
    disp('There is an error in Tc');
    Pi=[];Ti=[];Flag=0;
else
    %Compute choices and forks
    Tf = [];Pc = [];        
    for i=1:numT
        if size(find(Post(:,i)),1)>1
            Tf = [Tf i];
        end
    end
    for j=1:numP
        if size(find(Pre(j,:)),2)>1
            Pc = [Pc j];
        end
    end
    sc1 = 0;sc2 = 0;sc3 = 0;
    %Compute influence
    [Pi,Ti,sc1] = SimHPN_Influ(Pre,Post,Tc);
    
    if sc1==1             %SC1 is fulfilled?
        if isempty(Pc)     %SC2 is fulfilled?
           sc2 = 1;
        else
           pc = [];
           for i=1:size(Pc,2)
               postpc = find(Pre(Pc(i),:))';
               if all(ismember(postpc,Tc))
                   pc = [pc 1];
               else
                   pc = [pc 0];
               end  
           end
           if all(pc)
                sc2 = 1;
           end
        end
        
        if sc2==1
            if isempty(Tf)
                sc3 = 1;
            else
            tf = [];
            for i=1:size(Tf,2)
               post_tf = find(Post(:,Tf(i)))';
               ppost_tf = [];
               for j=1:size(post_tf,2)
                   ppost_tj = find(Pre(post_tf(j),:));
                   ppost_tf = [ppost_tf ppost_tj];
               end
               ppost_tf = unique(ppost_tf);
               if size(setdiff(ppost_tf,Tc),2)<=1
                   tf = [tf 1];
               else
                   tf = [tf 0];
               end  
            end
            if all(tf)
                   sc3 = 1;
            end
            end
            if sc3==1
                Flag=1;
            end
        end    
    end
%if Flag==1
%        disp('The timed net is NRC')
%elseif sc1==0
%        disp('The net is not NRC')
%else
%        disp('It is not possible to determine if the net is NRC')
%end
%    
end

    
  
return

