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

% SimHPN_control_onoff() standard ON-OFF controller for Choice-free TCPN

%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre, Post: Pre- and post- incidence matrix
% m0:        Initial marking
% lambda:    Transitions firing rates
% delta:     Sampling period
% mf:        Final marking
% fv:        The firing count vector, if fv == [], using the minimal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Output parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% s:  number of steps for reaching the final marking
% m:  |P|*(s+1) matrix, marking trajectory of places, m(:, k) is the
%     state in k_th step, m(:, 1) = m0, m(:, s+1) = mf
% u:  |T|*(s) matrix, control actions
% w:  |T|*(s) matrix, the controlled flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ s, m, u, w ] = SimHPN_control_onoff( Pre, Post, m0, lambda, delta, mf, fv )

% net structure verification 
% ON-OFF controller is only applicable for stucturally persistent TCPN if:
% (A1) The flow matrix C has full rank, OR
% (A2) The net is strongly connected and consistent.

% Norm2 distance threshold to final state
qdis = 1e-6;

[PlaceNum, TransNum] = size(Pre);
verified = 1;

% incidence matrix
C = Post - Pre;

% structures flags
spFlag = 1; 
frFlag = 1;
sconcFlag = 1;
consistentFlag = 1;

% check if the system is structurally persistent
for i = 1:PlaceNum
    tempRi = (Pre(i, :)> 0);
    if sum(tempRi) > 1
        spFlag = 0;
        fprintf('The net is not stucturally persistent, ON-OFF controller may be not minimal time here\n')
        break;
    end
end

if spFlag == 1
    fprintf('The net is stucturally persistent\n');
end

    
if rank(C) == min(size(C))
    fprintf('The net has full rank flow matirx\n');
else
    fprintf('The net does not have full rank flow matirx\n');
    frFlag = 0;
end

% check if the system is strongly connected and consistent
if (SimHPN_control_strConn(Pre, Post) == 1) 
    fprintf('The net is strongly connected \n');
    if (SimHPN_control_consChk(C) == 1)
        fprintf('The net is consistent \n');
    else
        consistentFlag = 0;
        fprintf('The net is not consistent \n');
    end
else
    fprintf('The net is not strongly connected \n');
    sconcFlag = 0;
end


if (size(fv, 1) == 0)
    % Compute minimal firing count vector
    minFireVec = SimHPN_control_minFv(C, m0, mf);    
    fprintf('\nCorresponding minimal firing count vector is: %s\n', mat2str(minFireVec,4));
else
    fprintf('\nUsing a given firing count vector:%s\n', mat2str(fv,4));
    minFireVec = fv;
end

if isequal(minFireVec, zeros(size(Pre,2), 1))
    s = 0;
    return;
end

% get the starting time tick
tStart = tic;

global tEscape;
% Compute and apply control actions

% the accumulative amount of firings of transitions
sumFv = zeros(TransNum, 1);

enab = zeros(TransNum, 1);
flow = zeros(TransNum, 1);
u = zeros(TransNum, 1);

flag = 1;  % if flag == 0, all the transitions have reached the total
          % required firing count
m(:, 1) = m0;
s = 1;

controlledFlow = [];
%stopInstant = zeros(TransNum,1);

stoppedFlag = zeros(TransNum,1);

while flag == 1
    flag = 0;
    for i = 1:TransNum       
        % compute enabling degree
        ip = find(Pre(:,i));
        ms = m(:, s);   % the marking in step s
        mi = ms(ip);   
        Prei = Pre(:, i);
        PreiNonZero = Prei(ip);
        wmi = mi./PreiNonZero;
        enab(i) = min(wmi);
        
        flow(i) = lambda(i) * enab(i);
        
        % compute the control law        
        if (sumFv(i) + flow(i) * delta > minFireVec(i)) && (sumFv(i) < minFireVec(i))   % The last firing of t_i
            u(i, s) = flow(i) - (minFireVec(i) - sumFv(i))/delta;
            flag = 1;
            %stopInstant(i) = stopInstant(i) + 1;
        elseif  sumFv(i) < minFireVec(i) % Totally ON
            u(i, s) = 0;
            flag = 1;
            %stopInstant(i) = stopInstant(i) + 1;
        else                             % Totally OFF 
            stoppedFlag(i) = 1;
            u(i, s) = flow(i, 1);             
        end
    end
    
    % if the state is closed enough to the final marking 
    dif = m(:, s) - mf;
    if dif' * dif < qdis
        flag = 0; % stop the process
    end

    if flag == 1
        % update the accumulative flow
        fc = (flow(:, 1) - u(:,s)) * delta;
        sumFv = sumFv + fc;

        controlledFlow(:, s) = flow(:, 1) - u(:,s);
        
        % update the marking
        s = s + 1;
        m(:, s) = m(:, s - 1) + C * fc;
    end
end

w = controlledFlow;

fprintf('stoppedFlag = %s\n', mat2str(stoppedFlag));

% nubmer of steps reaching final states
s = s -1;

% get ending time tick
tEscape = toc(tStart);
fprintf('Time consumed: %f millisecond, step: %d\n', tEscape * 1000, s);



