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

%THIS ALGORITHM COMPUTES A CONTROL LAW (Fm+g) FOR A TCPN ASSUMING ALL THE
%TRANSITIONS AS CONTROLLABLE, CONSIDERING ALMOST ALL CLASS(m0)  (alpha=1)
function [F,g]=SimHPN_control_All_T_contF2(C,L,mq,m0,alpha,wb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     COMPUTE BASIC PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
By=(null(C'))';
d=rank(C);
[np,nt]=size(C);

epsilon=0.001;
lbw=wb(1);
ubw=wb(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     COMPUTE VERTICES OF CLASS(M0) AND VERTICES RELATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building equalities and inequalities that characterize Class
Aineq=[-eye(np)];
bineq=[zeros(np,1)];  
Aeq=[By];
beq=[By*m0];
rAeq=size(Aeq,1);
%%%%% Compute all the possible crossing hyperplanes, i.e., probable
%%%%% vertices
Tri=combnk([1:np],np-rAeq);
[rTri,cTri]=size(Tri);

V=[];
R=[];
for i=1:rTri
    sA=[Aineq(Tri(i,:),:);Aeq];
    sb=[bineq(Tri(i,:));beq];
    if rank(sA)==np
        vi=sA\sb;
        if (isfinite(vi)==ones(np,1))&(sum(Aineq*vi<=bineq)==np)
            [rv,cv]=size(V);
            if cv==0
                V=[V,vi];
                R=[R,Tri(i,:)'];
            else
                eq=0;
                for j=1:cv
                    if V(1:np,j)==vi(1:np)
                        eq=1;
                    end;
                end;
                if eq==0
                    V=[V,vi];
                    R=[R,Tri(i,:)'];
                end;
           end;
        end;
    end;
end;
[rv,nv]=size(V);

%%% Reordering the vectors, so the first d+1 are linearly independent
indvect=[1];
depvect=[];
for i=2:nv
    if rank(V(:,[indvect,i]))>length(indvect)
        indvect=[indvect,i];
    else
        depvect=[depvect,i];
    end;
end;
V=V(:,[indvect,depvect]);
R=R(:,[indvect,depvect]);
        
%Matrices for the continuous dynamical system
A=zeros(np,np);
B=C;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     COMPUTATION  OF INPUTS AT VERTICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%First, compute the vectors orthogonal to facets
normalF=[];
nfacet=np;
mN=zeros(1,np,nv);
nofas=zeros(nv,1);
for i=1:nfacet
    Indexi=[];
    [I,Indexi]=find(R==i);
    Indexi=unique(Indexi);
    %First case, the facet has more than one vertex
    if (isempty(Indexi)==0)&(length(Indexi)>1)
        Res=[];
        for k=1:length(Indexi)-1
            Res=[Res,V(:,Indexi(1))];
        end;
        Facet=V(:,Indexi(2:length(Indexi)))-Res;
        %Compute the normal and save it
        ortho=null([Facet';By]);
        ortho=-ortho*sign(ortho'*(mq-V(:,Indexi(1))));
        normalF=[normalF,ortho];
        %Save the normal in the corresponding vertex matrix
        for k=1:length(Indexi)
            nofas(Indexi(k))=nofas(Indexi(k))+1;
            mN(nofas(Indexi(k)),1:np,Indexi(k))=(normalF(:,i))';
        end;
    %Second case, the facet has only one vertex
    elseif (isempty(Indexi)==0)&(length(Indexi)==1)
        ortho=null(By);
        ortho=-ortho*sign(ortho'*(mq-V(:,Indexi(1))));
        normalF=[normalF,ortho];        
        nofas(Indexi(1))=nofas(Indexi(1))+1;
        mN(nofas(Indexi(1)),1:np,Indexi(1))=(normalF(:,i))';
   %Third case, the facet is unfeasible
    else
        normalF=[normalF,zeros(np,1)];
    end;
end;

%%%%%%Secondly, compute the matrices for the inequality condition
[ns,nk]=size(mN(:,:,1)*B);
bineq=[mN(:,:,1)*A*V(:,1)];
Aineq=[-mN(:,:,1)*B];
for i=2:nv
    bineq=[bineq;mN(:,:,i)*A*V(:,i)];
    Aineq=[Aineq,zeros((i-1)*ns,nk);zeros(ns,nk*(i-1)),-mN(:,:,i)*B];
end;
bineq=-bineq;
Aineq=-Aineq;

% %%%%%% Add to the inequality, the condition: vertices are not equilibrium points
% Aneqp=[];
% bneqp=[];
% for i=1:nv
%     vA=zeros(1,nv*nk);
%     vA((i-1)*nk+1:i*nk)=ones(1,d)*mN(:,:,i)*B;
%     Aneqp=[Aneqp;vA];
%     vB=-epsilon-ones(1,d)*mN(:,:,i)*A*V(:,i);
%     bneqp=[bneqp;vB];
% end;
% Aineq=[Aineq;Aneqp];
% bineq=[bineq;bneqp];

%%%%%%Next, compute the matrices for the equality condition
%%%%%%Assuming that the first d+1 vertices are linearly independent
if nv>d+1
    mP=[V(:,d+2:nv);ones(1,nv-d-1)]'/[V(:,1:d+1);ones(1,d+1)]';
    Mu=[mP,-eye(nv-d-1)];
    
    if nk==1
        Aeq=Mu;
        beq=zeros(nv-d-1,1)
    else
        Aeq=[];
        for j=1:nv    
             subAeq=Mu(:,j);
             for i=2:nk
                 subAeq=[subAeq,zeros((i-1)*(nv-d-1),1);zeros((nv-d-1),i-1),Mu(:,j)];
             end;
             Aeq=[Aeq,subAeq];            
        end;
        beq=zeros((nv-d-1)*nk,1);      
    end;
else
    Aeq=[];
    beq=[];
end;



%%%%%%%%%%% Constraint for equilibrium point
 gama=V(:,1:(d+1))\mq;
 gama=[gama;zeros(nv-d-1,1)];
 Aeqg=[];
 for j=1:nv
     subAeqg=gama(j);
     for i=2:nk
         subAeqg=[subAeqg,zeros((i-1),1);zeros(1,i-1),gama(j)];
     end;
     Aeqg=[Aeqg,subAeqg];
 end;
 Aeqg=B*Aeqg;
 beqg=zeros(np,1);
 Aeq=[Aeq;Aeqg];
 beq=[beq;beqg];

% %Computing the upperbound for the input
[np,nt]=size(C);
Pre=-C.*(C<0);
Pre=Pre';
enab=zeros(nt,nv);
upb=[];
for i=1:nv
    for j=1:nt
        vj=find(Pre(j,:));Auxfj=Pre(j,:).*(V(:,i)*alpha+mq*(1-alpha))';
        enab(j,i)=min(Auxfj(vj));
    end;
    upb=[upb;diag(L)*enab(:,i)];
end;

% %Improving optimization (a larger flow upperbound is considered for vertices closer to m1 and m0)
% vdistance=[];
% for i=1:nv
%     vdistance=[vdistance,(V(:,i)-m0)'*(V(:,i)-m0)+(V(:,i)-mq)'*(V(:,i)-mq)];
%     [Y,I]=sort(vdistance,'descend');
% end;
% Jvector=[];
% for i=1:nv
%     Jvector=[Jvector,ones(1,nt)*find(i==I)*((ubw-lbw)/(nv-1))];
% end;
% Jvector=Jvector+ones(1,nt*nv)*lbw;

%Improving optimization (maximize the projection of the controlled flow on the direction to mf)
Ldir=[];
for i=1:nv
    Ldir=[Ldir -(mq-V(:,i))'*C];
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     SOLVING THE LPP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 options=optimset('Display','off');
 %option 1, the best in second example datosLevid2
%[U,FVAL,EXITFLAG]=fmincon(@(x)((x-upb)'*(x-upb)),zeros(nv*nk,1),Aineq,bineq,Aeq,beq,zeros(nv*nk,1),Jvector'.*upb,[],options);
 %option 2
%[U,FVAL,EXITFLAG]=quadprog(eye(nt*nv),zeros(nt*nv,1),Aineq,bineq,Aeq,beq,zeros(nv*nt,1),Jvector'.*upb,zeros(nv*nt,1),options);
%option 3
[U,FVAL,EXITFLAG]=linprog(Ldir',Aineq, bineq, Aeq, beq,zeros(nv*nt,1),10*ones(nv*nt,1), zeros(nv*nt,1),options);
%[U,FVAL,EXITFLAG]=fmincon(@(x)(Ldir*x),zeros(nv*nk,1),Aineq,bineq,Aeq,beq,zeros(nv*nk,1),Jvector'.*upb,[],options);

if EXITFLAG<0
     disp('NO SOLUTION FOUND');
     FVAL
 elseif EXITFLAG==0
     disp('MAXIMUM NUMBER OF ITERATIONS REACHED');
     FVAL
 end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     COMPUTE THE MATRIX GAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if EXITFLAG>0
    disp('CONTROL INPUTS SUCESSFULLY COMPUTED');
    FVAL
    Um=[];
    for i=1:nv
        Um=[Um,U((i-1)*nk+1:i*nk)];
    end;
    Sol=Um(:,1:d+1)/[V(:,1:d+1);ones(1,d+1)];
    F=Sol(:,1:np)
    g=Sol(:,np+1)
    if rank([A+B*F;By])==np
        disp('Equilibrium point is unique');
    else
        disp('Equilibrium point is not unique');
    end;
end;
