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

function f = SimHPN_resf(m, pre, post, la, eabs)
global nps nts;

% Computes the flow through transitions applying
% finite servers semantics (strong and weak enabling).

if m(2)<-.5
    3;
end

nps=size(pre,1);   % Number of places
nts=size(pre,2);	 % Number of transitions


% Compute the empty input places of each transition

[nemps,emps]=com_emps(m, pre, eabs); % nemps(i) stores the number of empty input places
% of transition i. emps is a matrix in which the row i stores the empty
% input places of transition t

% Compute the input flow to empty places

ifemps=com_ifemps(m, pre, post, eabs); % ifemps(i) is the input flow of place i

% Compute the flow through each transition

syms eqs; % symbolic flow equations

for i=1:nts
    if emps(i,1)==0  % No empty imput places. Strong enabling.
        eqs(i)=strcat('f',num2str(i),'=',num2str(la(i)));
    else				  % Empty places. Weak enabling.
        equ=strcat('f',num2str(i),'=min(',num2str(la(i)));
        for j=1:nemps(i) % Minimum flow of input places
            curif=strcat('(',char(ifemps(emps(i,j))),')');
            % if the net is join free
            %     curpre=num2str(pre(emps(i,j),i)); can be used
            % else a flow splitting among output transitions
            % of the empty place is done
            curpre=num2str(sum(pre(emps(i,j),:)));
            curter=strcat(curif,'/',curpre);
            equ=strcat(equ,',min(',curter);
        end
        for j=1:nemps(i)+1 % Close brackets
            equ=strcat(equ,')');
        end
        eqs(i)=equ;
    end
end

% Flow equations system

feq=strcat(char(eqs(1)));

for i=2:nts
    feq=strcat(feq,',',char(eqs(i)));
end

% Flow values

if nts==1
    solf.f1=solve(feq);
else
    try
        solf=solve(feq);
    catch
        disp('Cannot compute immediate throughput. Check loops of empty places!');
        for i=1:nts
            eval(strcat('solf.f',num2str(i),'=0;'));
        end
    end
end


% Extract numeric values

if size(solf.f1,1)>1  % Multiple solutions for weak enabling due to empty marking
    for i=1:nts
        curf=eval(strcat('solf.f',num2str(i)));
        % If there is a null throughput take it, else take any of the solutions
        valf=1;j=1;
        while valf~=0 && j<=size(curf,1)
            valf=eval(strcat('curf(',num2str(j),')'));
            j=j+1;
        end
        if valf==0
            f(i)=0;
        else
            f(i)=eval(strcat('curf(',num2str(1),')'));
        end
    end
else
    for i=1:nts
        f(i)=eval(strcat('solf.f',num2str(i)));
    end
end

% No more symbolic variables

try
    f=eval(f)';
catch		% f is already non-symbolic
    f=f';
end


%%%%% Compute empty places function %%%%%

function [nemps,emps] = com_emps(m, pre, eabs)
global nps nts ;

nemps=zeros(nts,1);
emps=zeros(nts,1);

for i=1:nps
    if m(i)<=eabs % Marking considered 0 if less than the allowed error
        for j=1:nts
            if pre(i,j)>0 % i is an empty place of j
                nemps(j)=nemps(j)+1;
                emps(j,nemps(j))=i;
            end
        end
    end
end

%%%%% Compute input flow of empty places function %%%%%

function ifemps = com_ifemps(m, pre, post, eabs)
global nps nts ;

syms ifemps

for i=1:nps
    ifemps(i)='0';
    if m(i)<=eabs % Marking considered 0 if less than the allowed error
        for j=1:nts
            if post(i,j)>0
                curpost=num2str(post(i,j));
                ifemps(i)=strcat(char(ifemps(i)),'+f',num2str(j),'*',curpost);
            end
        end
    end
end
