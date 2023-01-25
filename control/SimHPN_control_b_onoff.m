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

% SimHPN_control_b_onoff(): balanced ON-OFF controller, it is assumed that m0 > 0

%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre, Post: Pre- and post- incidence matrix
% m0:        Initial marking
% lambda:    Transitions firing rates
% delta:     Sampling period
% mf:        Final marking
% fv:        The firing count vector, if fv == [], using a minimal one
% d_threshold: used to distinguish fast and slow transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Output parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% s:  number of steps for reaching the final marking
% m:  |P|*(s+1) matrix, marking trajectory of places, m(:, k) is the
%     state in k_th step, m(:, 1) = m0, m(:, s+1) = mf
% u:  |T|*(s) matrix, control actions
% w:  |T|*(s) matrix, the controlled flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ s, m, u, w ] = SimHPN_control_b_onoff( Pre, Post, m0, lambda, delta, mf, fv, d_threshold )

[placeNum, transNum] = size(Pre);
C = Post - Pre;

% Norm2 distance threshold to final state
qdis = 1e-6;

if (size(fv, 1) == 0)
    % Compute minimal firing count vector
    fv = SimHPN_control_minFv(C, m0, mf);    
    fprintf('\nCorresponding minimal firing count vector is: %s\n ', mat2str(fv,4));  
else
    fprintf('\nUsing a given firing count vector: %s\n', mat2str(fv,4));
end


% get the starting time tick
tStart = tic;

fv_remain = fv;
s = 0; % the number of time step
mcur = m0;
m = m0;
u = [];
w = [];

% the longest path for each transition to get tokens, for extension
longestpath = ones(transNum,1);

% the number of time steps in which D is not decreased
non_decrease = zeros(transNum,1);

% estimated number of stpes
estimated_step = [];

%initial D value
D_cur = zeros(transNum,1);

% get the confliting matrix
cm = SimHPN_control_cflBalance(Pre);

t_enab = ones(transNum,1);
%initiate t_eabl: if t_j is not in a conflit, t_enab[j] = 1, else 0
for i = 1 : size(cm,1)
    for j = 1 : size(cm,2)
        if cm(i,j) ~= 0
            t_enab(j) = 0;
        end
    end
end

% check the estimated number of stpes for each group of conflicting
% transitions in the initial state
    
for i = 1 : size(cm,1)
    for j = 1 : size(cm,2)
        if cm(i,j) ~= 0
            f = lambda(j) * SimHPN_control_enabling(Pre, j, m0) * delta;
            if f == 0
                estimated_step(i, j) = inf;
            else
                estimated_step(i, j) = fv(j) / f;
            end
        else
            estimated_step(i,j) = 0;
        end
        
    end
end

% divide each group of transitions into fast and slow
for (i = 1 : size(estimated_step,1))        
    %set: t_enab[j] = 1 if t_j is a fast transition, else 0         
    row = estimated_step(i, :);

    while ~isempty(row(row > 0))
        minstep = min(row(row>0)); % minimal number of step of the row
        temp = find(row == minstep);
        tj = temp(1); % get the transition which has the minimal step

        t_enab(tj) = 1; % tj should be enabled
        row(tj) = 0; 
        
        nextminstep = min(row(row > 0));
        if ceil(minstep) * d_threshold < ceil(nextminstep)
            break;
        end
    end
end

for j = 1 : transNum
    if t_enab(j) == 0
        D_cur(j) = inf;  %%% ???
        fv_remain(j) = 0;
    end
end

while ((mcur - mf)' * (mcur - mf) > qdis)   
    % apply the control
    [mcur, ucur, wcur] = SimHPN_control_b_onoff_LPP(Pre, Post, mcur, lambda, delta, fv_remain, cm);
    m = [m, mcur];
    w = [w, wcur];       
    u = [];  % computed in the calling function
    fv_remain = fv_remain - wcur * delta;

    % update the fast and slows sets of each group of conflicting
    % transitions
    
    % get the flow of the next step !!!
    [mtemp, utemp, wtemp] = SimHPN_control_onoff_plus_LPP( Pre, Post, mcur, lambda, delta, fv_remain); % may use a more specific function to save time
    
    estimated_step = [];
    
    % update the estimated number of steps
    for i = 1 : size(cm,1)
        for j = 1 : size(cm,2)
            if cm(i,j) ~= 0
                % if t_j is a fast transition
                if t_enab(j) == 1
                    f = wtemp(j) * delta;
                    if f == 0
                        estimated_step(i, j) = inf;
                    else
                        estimated_step(i, j) = fv_remain(j) / f;
                    end
                else % if t_j is a slow transition
                    f = lambda(j) * SimHPN_control_enabling(Pre, j, mcur) * delta;
                    if f == 0
                        estimated_step(i, j) = inf;
                    else
                        estimated_step(i, j) = fv(j) / f;
                    end
                end
            else
                estimated_step(i, j) = 0;
            end
        end
    end
    
    %estimated_step 

    for i = 1 : size(cm,1)
        for j = 1 : size(cm, 2)
            if t_enab(j) == 1
                h = estimated_step(j);
                break;
            end
        end
        for j = 1 : size(cm,2) 
            if t_enab(j) == 0
                s_j = estimated_step(j);              

                D_last =  D_cur(j);
                D_cur(j) = s_j / h;
                if D_cur(j) >= D_last
                    non_decrease(j) = non_decrease(j) + 1;                    
                end
                if D_cur(j) <= d_threshold || (D_cur(j) >= D_last && non_decrease(j) >= longestpath(j))                                        
                    t_enab(j) = 1; 
                    fv_remain(j) = fv(j);
                end
            end
        end
    end

    % number of time steps
    s = s + 1;
    
end

global tEscape;
% get ending time tick
tEscape = toc(tStart);
fprintf('Time consumed: %f millisecond, step: %d\n', tEscape * 1000, size(m,2) - 1);
