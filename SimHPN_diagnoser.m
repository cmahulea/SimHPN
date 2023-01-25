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

function [BM,BMV,delta] = SimHPN_diagnoser(Pre,Post,m0,obsT,obsW,obsA,Tf)
% Pre and Post - incidence matrices
% m0 - the initial marking
% obsT - number of observable transitions (first obsT transitions will be observable and the other one not)
% obsW - the observed sequence of transitions
% obsA - the observed cuantities of firing correspopnding to each
% transition in obsW

if isempty(obsT)
    errordlg('Please introduce the number of observable transitions.');
    return;
end
if isempty(obsW)
    errordlg('Please introduce the observed word.');
    return;
end
if isempty(obsA)
    errordlg('Please introduce for each transition of the observed word its firing quantity.');
    return;
end
if (length(obsW) ~= length(obsA))
    errordlg('Number of observed transitions should be equal with the observed quantity');
    return
end

Cin = Post - Pre;
Tu = setdiff(1:size(Cin,2),1:obsT);
Cu = Cin(:,Tu);
%%%%%%%%%%%%enabling bounds of transitions%%%%%%%%%%%%%%
bounds = zeros(length(Tu),1);
for i = 1 : length(Tu)
    IN.obj = [zeros(1,size(Cin,2)) -1];
    IN.A = [-Cin Pre(:,Tu(i))];
    IN.B = m0;
    IN.A = [IN.A ; -eye(size(IN.A,2))];
    IN.B = [IN.B ; zeros(size(IN.A,2),1)];
    OUT=cddmex('solve_lp',IN);
    if (OUT.how == 1)
        bounds(i) = OUT.xopt(length(OUT.xopt));
        fprintf(1,'Enabling bound of t%d is %d\n',Tu(i),bounds(i));
    else
        error('Error when computing the structural enabling bounds of the unobservable transitions. Maybe the net is not bounded!');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = [eye(length(m0)) -Cu];
B = m0;
A = [A; -eye(length(m0)) Cu];
B = [B; -m0];
A = [A; -eye(length(m0)) zeros(length(m0),length(Tu))];
B = [B;zeros(length(m0),1)];
A = [A; zeros(length(Tu),length(m0)) -eye(length(Tu))];
B = [B;zeros(length(Tu),1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [A;zeros(length(bounds),length(m0)) eye(length(bounds))];
B = [B;bounds];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=struct('A',A,'B',B);
temp =  cddmex('extreme',H);

BMV{1} = abs(round(temp.V*10000)/10000);%intialiy, the set of consistent markings
BM{1}= H;

fprintf(1,'\n***********************************************************\n');
fprintf(1,'For empty word:\n');
fprintf(1,'***********************************************************\n\n');

fprintf(1,'Vertices of the set of consistent markings (vertices of \\bar Y(m_0,w))\n');
temp = BMV{1};
for k = 1 : size(temp,1)
    fprintf(1,'e_%d = %s''\n',k,mat2str(temp(k,:)));
end

fprintf(1,'\nComputing the diagnosis states\n');

for k = 1 : length(Tf)
    index = zeros(1,length(m0) + length(Tu));
    faults = Tf{k};
    for ll = 1 : length(faults)
        index(length(m0)+faults(ll)-obsT)=1;
    end
    IN = BM{1};
    IN.obj = index;
    OUT=cddmex('solve_lp',IN);
    if (OUT.how == 1)
        li = OUT.objlp;
    elseif (OUT.how == 6)
        li = Inf;
    else
        error('Unable to solve LPP to compute li!!')
    end
    IN.obj = -index;
    OUT=cddmex('solve_lp',IN);
    if (OUT.how == 1)
        ui = -OUT.objlp;
    elseif (OUT.how == 6)
        ui = Inf;
    else
        error('Unable to solve LPP to compute ui!!')
    end
    if (ui == 0)
        delta{1}.k = 'N';
    elseif ((li == 0) && (ui > 0))
        delta{1}.k = 'U';
    elseif (li > 0)
        delta{1}.k = 'F';
    end
    fprintf(1,'Diagnosis state for Tf^{%d}: %s\n',k,delta{1}.k);
end

for i = 1 : length(obsW)    
    fprintf(1,'\n***********************************************************\n');
    str = 'Observed sequence: ';
    for k = 1 : i
        str = sprintf('%s t%d(%4.2f)',str,obsW(k),obsA(k));
    end
    disp(str);

    H = BM{i};
    H.A = [H.A; -eye(length(m0)) zeros(length(m0),length(Tu))];
    H.B = [H.B; -obsA(i) * Pre(:,obsW(i))];
    tic;
    temp = cddmex('extreme',H);
    time_c = toc;
    fprintf(1,'***********************************************************\n\n');
    fprintf(1,'Vertices of the set of consistent markings after cutting\n');
    temp = abs(round(temp.V*10000)/10000);    

    for k = 1 : size(temp,1)
        fprintf(1,'e_%d = %s''\n',k,mat2str(temp(k,:)));
    end
    vertices = temp;
    E = [];
    tic;
    for j = 1 : size(vertices,1)
        current = vertices(j,:);
        % variables are: [m \sigma_u]
        % first constraints: m = \tilde{m} + Cu \sigma_u + \alpha C(:,t)
        A = [eye(length(m0)) -Cu];
        B = [current(1:length(m0))'+obsA(i)*Cin(:,obsW(i))];
        A = [A; -eye(length(m0)) Cu];
        B = [B ; -current(1:length(m0))'-obsA(i)*Cin(:,obsW(i))];
        % second constraint: m, \sigma_u, \tilde{m}, \sigma_u' \geq 0
        temp = size(A,2);
        A = [A; -eye(temp)];
        B = [B; zeros(temp,1)];
        %%%%%%%%%enabling bounds%%%%%%
        A = [A;zeros(length(bounds),length(m0)) eye(length(bounds))];
        B = [B;bounds];
        %%%%%%
        H=struct('A',A,'B',B);
        temp =  cddmex('extreme',H);
        temp = abs(round(temp.V(:,1:length(m0)+length(Tu))*10000)/10000);
        for k = 1 : size(temp,1)
            temp(k,length(m0)+1:end) = temp(k,length(m0)+1:end)+ current(length(m0)+1:end);
        end
        E = [E ; temp];
    end
    V = struct('V',E);
    next1 = cddmex('hull',V);
    time_c = time_c + toc;
    next2 = cddmex('extreme',next1);
    BMV{i+1}= abs(round(next2.V*10000)/10000);
    BM{i+1} = next1;
    fprintf(1,'***********************************************************\n');
    fprintf('Vertices of the set of consistent markings (vertices of \\bar Y(m_0,w))\n');
    temp = BMV{i+1};
    for k = 1 : size(temp,1)
        fprintf(1,'e_%d = %s''\n',k,mat2str(temp(k,:)));
    end
    fprintf(1,'Computation time: %20.2f seconds',time_c);
    fprintf(1,'\nComputing the diagnosis states\n');
        for k = 1 : length(Tf)
        index = zeros(1,length(m0) + length(Tu));
        faults = Tf{k};
        for ll = 1 : length(faults)
            index(length(m0)+faults(ll)-obsT)=1;
        end
        IN = BM{i+1};
        IN.obj = index;
        OUT=cddmex('solve_lp',IN);
        if (OUT.how == 1)
            li = OUT.objlp;
        elseif (OUT.how == 6)
            li = Inf;
        else
            error('Unable to solve LPP to compute li!!')
        end
        IN.obj = -index;
        OUT=cddmex('solve_lp',IN);
        if (OUT.how == 1)
            ui = -OUT.objlp;
        elseif (OUT.how == 6)
            ui = Inf;
        else
            error('Unable to solve LPP to compute ui!!')
        end
        if (ui <= eps)
            delta{i+1}.k = 'N';
        elseif ((li <= eps) && (ui > eps))
            delta{i+1}.k = 'U';
        elseif (li > eps)
            delta{i+1}.k = 'F';
        end
        fprintf(1,'Diagnosis state for Tf^{%d}: %s\n',k,delta{i+1}.k);
    end
end
return
