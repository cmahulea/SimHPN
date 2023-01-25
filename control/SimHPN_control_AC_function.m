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

function [timeline,marking,fflow,caction,F,g]=SimHPN_control_AC_function(Pre,Post,L,m0,mq,delta,ftime,wb);

global tEscape;
tStart = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        COMPUTATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=Post-Pre;PreM=Pre';
[np,nt]=size(C);
%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROL LAW SYNTHESIS
alphai=0.9;
if exist('wb','var')==0
    wb=[1;1];
else
    if (size(wb,1)~=2)|(size(wb,2)~=1)
        wb=[1;1];
    end;
end;
[F,g]=SimHPN_control_All_T_contF2(C,L,mq,m0,alphai,wb);
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZING
m=m0;
timeline=[];
marking=[];
fflow=[];
caction=[];
disp('SIMULATING CLOSED-LOOP SYSTEM...');
%for k=1:floor(ftime/delta)
 k=0;
 while (m-mq)'*(m-mq)>=1e-6
     k=k+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% FLOW AND ENABLING DEGREE COMPUTATION
    enab=zeros(nt,1);
    for j=1:nt
        vj=find(PreM(j,:));Auxfj=PreM(j,:).*m';
        enab(j)=min(Auxfj(vj));
    end;
    flow=diag(L)*enab;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONTROL LAW  EVALUATION
%    alpha=1-min(m./mq);
%     if alpha>0.95
%         w=flow;
%     else
         w=F*m+g;
%    end;
    Index=find(w>0);vectcomp=flow(Index)./w(Index);beta=min(vectcomp);
    input=flow-beta*w;
    %%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATE UPDATING
    m=m+C*(flow-input)*delta;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% SAVING DATA
    timeline=[timeline, delta*k];
    marking=[marking,m];
    fflow=[fflow,flow-input];
    caction=[caction,input];   
end;

% get ending time tick
tEscape = toc(tStart);
fprintf('Time consumed: %f millisecond\n', tEscape * 1000);

disp('SIMULATION DONE');

