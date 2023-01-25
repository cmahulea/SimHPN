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

function [MeanMarking,MeanThroughput,Tline]=SimHPN_SimHyb(demo,Pre,Post,L,tipo,m0,nsim,tsim,seman,modo,delta)


if (nargin < 1)
    errordlg('This is a demo version. Please contact us by email at simhpn@unizar.es for the full version.','SimHPN Toolbox for MATLAB');
    return;
end
if (demo ~= -21)
    errordlg('This is a demo version. Please contact us by email at simhpn@unizar.es for the full version.','SimHPN Toolbox for MATLAB');
    return;
end
if (nargin ~= 11)
    errordlg('The number of input parameters is not correct!.','SimHPN Toolbox for MATLAB');
    return;
end

if (seman == 2)
    seman = 3;
end

Tline=[];MeanMarking=[];MeanThroughput=[];
%Tiempos=[];Datos=[];
%Checking for errors
errors=1;
if (size(Post,1)~=size(Pre,1))||(size(Post,2)~=size(Pre,2))
    errordlg('Matrices Post and Pre do not have the same dimensions','Data error');
    return;
elseif (size(L,2)~=1)||(size(L,1)~=size(Pre,2))
    errordlg('Vector lambda does not have the correct dimension','Data error');
    return;
elseif (size(tipo,2)~=1)||(size(tipo,1)~=size(Pre,2))
    errordlg('Vector type does not have the correct dimension','Data error');
    return;
elseif (size(m0,2)~=1)||(size(m0,1)~=size(Pre,1))
    errordlg('Vector m0 does not have the correct dimension','Data error');
    return;
elseif min(min(Pre))<0
    errordlg('Matrix Pre has negative values','Data error');
    return;
elseif min(min(Post))<0
    errordlg('Matrix Post has negative values','Data error');
    return;
elseif min(m0)<0
    errordlg('Vector m0 has negative values','Data error');
    return;
elseif length(find((tipo=='q')|(tipo=='d')|(tipo=='c')))<length(tipo)
    errordlg('Type of some transitions wrongly defined','Data error');
    return;
elseif min(L)<=0
    errordlg('Vector lambda has non-positive values','Data error');
    return;
elseif min(max(Pre',[],2))==0
    errordlg('According to Matrix Pre, some transition does not have input places','Data error');
    return;
end;

if (length(find((tipo(:)=='q')|(tipo(:)=='d')))<1)&&(nsim>1)
    nsimt=1;
else
    nsimt=nsim;
end;

%%%%%%%%%%%%%% PROCEED WITH THE SIMULATION %%%%%%%%%%%%%%
C = Post - Pre;
modo=modo-1;
[np,nt]=size(C);
PreM=Pre';PreM(PreM<=0)=-1;PreM=ones(size(PreM))./PreM;PreM(PreM<0)=0;
PM=[1:1:np];
TM=[1:1:nt];
%Relabeling transitions for product semantics
if seman==3
    tipo(find(tipo=='c'))='p';
    tipo(find(tipo=='d'))='m';
end;

%Computing sets of discrete and continuous transitions
sdt=find((tipo(:)=='q')|(tipo(:)=='d')|(tipo(:)=='m'));ndt=length(sdt);
sct=find((tipo(:)=='c')|(tipo(:)=='p'));nct=length(sct);
%Setting initial samplings exist('delta')
if (delta>0)&&(modo==0)
    dt=delta;
    dmuestreo=dt;
else
    dt=tsim;
    dmuestreo=tsim;
end;
%%%%%%%%%%%%%%%%%%% STARTING SIMULATIONS %%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Simulation progress','Name','SimHPN Toolbox');
ph = findobj(h, 'type', 'patch');
set(ph, 'FaceColor', [0 0 1], 'EdgeColor', [0 0 1]);

for l=0:nsimt
    %Simulation zero is only for unknown-fixed-delay operation mode exist('delta')
    if (l==0)&&(modo==1||modo==2||delta>0)
        l=1;
    end;
    
    
    %Initializing variables
    md=m0(1:np);count=zeros(nt,1);enab=zeros(nt,1);
    reloj=10*tsim*ones(nt,1);
    for i=1:nt
        enab(i)=0;
        vi=find(PreM(i,:));Auxfi=PreM(i,:).*md';
        if (tipo(i)=='q')||(tipo(i)=='d')||(tipo(i)=='c')
            enab(i)=min(Auxfi(vi));
        elseif (tipo(i)=='p')
            enab(i)=prod(Auxfi(vi));
        elseif (tipo(i)=='m')
            enab(i)=prod(floor(Auxfi(vi)));
        end;
        if (enab(i)<1)||(tipo(i)=='c')||(tipo(i)=='p')
            reloj(i)=10*tsim;
        elseif (enab(i)>=1)&&((tipo(i)=='q')||(tipo(i)=='d')||(tipo(i)=='m'))
            %Transition enabled, asign value for the reloj
            if (tipo(i)=='q')
                reloj(i)=1/L(i);
            elseif (tipo(i)=='d')||(tipo(i)=='m')
                reloj(i)=exprnd(1/(L(i)*floor(enab(i))));
            end;
        end;
        if (tipo(i)=='c')||(tipo(i)=='p')
            count(i)=L(i)*enab(i);
        end;
    end;
    time=0;k2=0;
    timeline=[time];Trayd=[md(PM);count(TM)];
    
    %%%%%%%%%%%%%%%%%%% SIMULATION CYCLE %%%%%%%%%%%%%%%%%%%
    while time<=tsim+dt
        %Compute the enabling degrees
        for i=1:nt
            enab(i)=0;
            vi=find(PreM(i,:));Auxfi=PreM(i,:).*md';
            if (tipo(i)=='q')||(tipo(i)=='d')||(tipo(i)=='c')
                enab(i)=min(Auxfi(vi));
            elseif (tipo(i)=='p')
                enab(i)=prod(Auxfi(vi));
            elseif (tipo(i)=='m')
                enab(i)=prod(floor(Auxfi(vi)));
            end;
        end;
        %Compute the next sampling
        if ((l==0)&&(modo==0))||(modo==1)||(modo==2)
            %if there exists an enabled transition
            if max(enab)>0
                dt=0.1/max(diag(L)*enab);
                dmuestreo=min(dt,dmuestreo);
                if (modo==1)||(modo==2)
                    dt=min([dt,min(reloj)]);
                end;
            else %if none transition is enabled, then there is a deadlock
                dt=tsim+dt-time;
                timeline=[timeline;tsim];
                Trayd=[Trayd [md(PM);count(TM)]];
            end;
        end;
        %Firing the continuous transitions
        if nct>0
            flujo=zeros(nt,1);
            for i=1:nct
                flujo(sct(i))=L(sct(i))*enab(sct(i))*dt;
                count(sct(i))=L(sct(i))*enab(sct(i));
            end;
            md=md+C*flujo;
        end;
        
        %Actualize relojs
        dif=reloj-ones(nt,1)*time;
        for i=1:ndt
            enab(sdt(i))=0;
            vi=find(PreM(sdt(i),:));Auxfi=PreM(sdt(i),:).*md';
            if (tipo(sdt(i))=='q')|(tipo(sdt(i))=='d')
                enab(sdt(i))=min(floor(Auxfi(vi)));
            elseif (tipo(sdt(i))=='m')
                enab(sdt(i))=prod(floor(Auxfi(vi)));
            end;
            if enab(sdt(i))<1
                %Transition not enabled = infinite reloj
                reloj(sdt(i))=time+10*tsim;
            elseif (enab(sdt(i))>=1)&(dif(sdt(i))>2*tsim)
                %Transition newly enabled, assign value to reloj
                if (tipo(sdt(i))=='q')
                    reloj(sdt(i))=time+1/L(sdt(i));
                elseif (tipo(sdt(i))=='d')|(tipo(sdt(i))=='m')
                    reloj(sdt(i))=time+exprnd(1/(L(sdt(i))*floor(enab(sdt(i)))));
                end;
            end;
        end;
        
        %Firing discrete transitions
        dif=reloj-ones(nt,1)*time;dif(sct)=10*tsim;
        if (min(dif)<=0)
            k=find(dif<=0);
            %Saving data if modo==1 or 2
            if (modo==1)|(modo==2)
                timeline=[timeline;time];
                Trayd=[Trayd [md(PM);count(TM)]];
            end;
            %Firing each discrete transition
            while length(k)>=1
                %Randomly chosing the transition to be fired
                if length(k)==1
                    j=1;
                else
                    j=ceil(rand*length(k));
                end;
                %If the enabling is >=1 and its time to be fired, then fire transition k(j)
                if (enab(k(j))>=1)%&(reloj(k(j))<=time)
                    mdold=md;
                    Parik=zeros(nt,1);Parik(k(j))=1;
                    md=md+C*Parik;
                    count(k(j))=count(k(j))+1;
                    
                    %Computing deterministic transitions in conflict
                    aux1=find(Pre(:,k(j))); %input places
                    aux2=find(sum(Pre(aux1,:),1)); %ouput transitions
                    aux3=aux2(find(tipo(aux2)=='q')); %deterministic transitions
                    %Actualize relojs
                    for i=1:ndt
                        vi=find(PreM(sdt(i),:));Auxfi=PreM(sdt(i),:).*md';
                        if (tipo(sdt(i))=='q')|(tipo(sdt(i))=='d')
                            enab(sdt(i))=min(floor(Auxfi(vi)));
                        elseif (tipo(sdt(i))=='m')
                            enab(sdt(i))=prod(floor(Auxfi(vi)));
                        end;
                        %Transition not enabled
                        if enab(sdt(i))<1
                            reloj(sdt(i))=time+10*tsim;
                            %Transition fired, newly enabled or deterministic in conflict
                        elseif (enab(sdt(i))>=1)&&((sdt(i)==k(j))||(dif(sdt(i))>2*tsim)||(~isempty(find(aux3==sdt(i)))))
                            if (tipo(sdt(i))=='q')
                                reloj(sdt(i))=time+1/L(sdt(i));
                            elseif (tipo(sdt(i))=='d')|(tipo(sdt(i))=='m')
                                reloj(sdt(i))=time+exprnd(1/(L(sdt(i))*enab(sdt(i))));
                            end;
                            %Transition remaining enabled
                        elseif (enab(sdt(i))>=1)&(dif(sdt(i))>0)&(dif(sdt(i))<2*tsim)&((tipo(sdt(i),:)=='d')|(tipo(sdt(i))=='m'))
                            enabold=0;
                            Auxfiold=PreM(sdt(i),:).*mdold';
                            if (tipo(sdt(i))=='d')
                                enabold=min(floor(Auxfiold(vi)));
                            elseif (tipo(sdt(i))=='m')
                                enabold=prod(floor(Auxfiold(vi)));
                            end;
                            reloj(sdt(i))=time+(reloj(sdt(i))-time)*enabold/floor(enab(sdt(i)));
                        end;
                    end;
                    %Eliminate from the list, transitions fired and
                    %deterministic transitions in conflict
                    k(j)=[];
                    if ~isempty(k)
                        for i=1:ndt
                            if (~isempty(find(aux3==sdt(i))))&&(~isempty(find(k==sdt(i))))
                                k(find(k==sdt(i)))=[];
                            end;
                        end;
                    end;
                end;
                %Eliminating from the list, transitions not enabled
                k(find(enab(k)<=0))=[];
            end;
            %Saving data if modo==1 or 2
            if (modo==1)||(modo==2)
                timeline=[timeline;time+dt/10];
                Trayd=[Trayd [md(PM);count(TM)]];
            end;
        end;
        %Actualize the time
        if ((modo==1)||(modo==2))&&(ndt==nt)
            dt=min(reloj)-time;
        end;
        time=time+dt;k2=k2+1;
        %Saving data every 10 sampling periods
        if (l>0)&&(mod(k2,10)==0)&&~((modo==1||modo==2)&&(ndt==nt))
            timeline=[timeline;time];
            Trayd=[Trayd [md(PM);count(TM)]];
        end;
        
    end;%%%%%%%%%%%%%%%%%%% END SIMULATION CYCLE %%%%%%%%%%%%%%%%%%%
    
    %Saving last simulation data
    if (l>0)
        timeline=[timeline;time];
        Trayd=[Trayd [md(PM);count(TM)]];
    end;
    
    %Assign the value for dt if modo==0 and compute the Sum of
    %trajectories
    [timeline indices J]=unique(timeline);
    Trayd=Trayd(:,indices);
    if modo==0
        if l==0
            dt=dmuestreo;
            fprintf(0,'Sampling settled as: %s seconds',num2str(dt));
        elseif l==1
            Sumad=Trayd';
        elseif l>1
            Sumad=Sumad+Trayd';
        end;
    elseif (modo==1)||(modo==2)
        if l==1
            tsim=min(tsim,timeline(length(timeline)));
            dt=dmuestreo;
            ddata=10*dt;
            basetime=ddata*[0:1:floor(tsim/ddata)];
            Sumad=zeros(length(basetime),size(Trayd,1));
        end;
        if (l>=1)&&length(timeline)>2
            for j=1:size(Trayd,1)
                Sumad(:,j)=Sumad(:,j)+interp1(timeline,Trayd(j,:),basetime');
            end;
        end;
    end;
    
    %Saving trajectories if modo==2
%     if (l==1)&(modo==2)
%         Datos=Trayd;
%         Tiempos=timeline';
%     elseif (l>1)&(modo==2)
%         [n1,n2,n3]=size(Datos);
%         nn=length(timeline);
%         if n2>nn
%             Tiempos=[Tiempos;[timeline',-1*ones(1,n2-nn)]];
%             Datos(:,:,l)=[Trayd,-1*ones(n1,n2-nn)];
%         else
%             Tiempos(:,n2+1:nn)=-1*ones(l-1,nn-n2);
%             Datos(:,n2+1:nn,:)=-1*ones(n1,nn-n2,l-1);
%             Tiempos=[Tiempos;timeline'];
%             Datos(:,:,l)=Trayd;
%         end;
%     end;
    %Simulation progress
    if ((mod(l,ceil(nsimt/40))==0)&&(l>0))||(l==nsimt)
        waitbar(l/nsimt)
    end;
    
end; %%%%%%%%%%%%%%%%%%% END SIMULATIONS %%%%%%%%%%%%%%%%%%%
if ~isempty(h)
    close(h)
end
%Set default values for some ouput variables
% if modo~=2
%     Tiempos=[];
%     Datos=[];
% else
%     assignin('base','C',C);
%     assignin('base','Pre',Pre);
%     assignin('base','lambda',L);
%     assignin('base','type',tipo);
%     assignin('base','semantics',seman);
%     assignin('base','M0',m0);
%     assignin('base','Simulation_parameter',[nsim;tsim;dt]);
%     assignin('base','Data',Datos);
%     assignin('base','Timelines',Tiempos);
% end;
if (modo==1||modo==2)
    timeline=basetime';
end;
%Computing average markings and throughputs
ntimew=10;
Tline=timeline;
MeanMarking=Sumad(:,1:length(PM))/nsimt;
MeanCount=Sumad(:,length(PM)+1:length(PM)+length(TM))/nsimt;
MeanThroughput=MeanCount;
sdtTM=find((tipo(TM)=='q')|(tipo(TM)=='d')|(tipo(TM)=='m'));
if ~isempty(sdtTM)
    if length(Tline)>ntimew
        MeanThroughput(:,sdtTM)=(MeanCount(:,sdtTM)-[zeros(ntimew,length(sdtTM));MeanCount(1:length(MeanCount(:,1))-ntimew,sdtTM)])/(ntimew*10*dmuestreo);
    else
        warndlg('Simulation time is too small for estimating throughputs','Simulation mode');
    end;
end;
