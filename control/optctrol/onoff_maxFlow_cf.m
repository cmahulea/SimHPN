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

% Drive the system to a equilbrium state with maximal flow, in minimal time
% WARNING: it is only for choice-free net systems at this moment 

%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre, Post: Pre- and post- incidence matrix
% m0:        Initial marking
% lambda:    Transitions firing rates
% delta:     Sampling period
% fv:        The firing count vector
% tj:        The transition to optimized
% maxFlow:   The maximal of flow of transitions tj in optimal steady states 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s, mTrajectory, uTrajectory, wTrajectory] = onoff_maxFlow_cf(Pre, Post, m0, lambda, delta, fv, tj, optFlow)

% the accumulative amount of firings of transitions
%sumFv = zeros(TransNum, 1);

[PlaceNum, TransNum] = size(Pre);

% computer the P,T-semiflow

[psemif, tsemif] = SimHPN_ptsemis(Post-Pre);

if size(tsemif,2) ~= 1
    fprintf('non CF net system\n');
else
    fprintf('the unique T-semfilow: %s\n', mat2str(tsemif,4));
end

intialFv = fv;

% remained fv that should be fired
reFv = fv;
sumFv = zeros(TransNum,1);

enab = zeros(TransNum, 1);
flow = zeros(TransNum, 1);
w = zeros(TransNum, 1);
u = zeros(TransNum, 1);

C = Post - Pre;

% the minimal makrings 
minMark = zeros(PlaceNum, 1);
for i = 1:PlaceNum
    j = find(Pre(i,:), 1); 
    minMark(i) = ((optFlow * tsemif(j) / tsemif(tj)) / lambda(j)) * Pre(i,j);
end

m = m0;
s = 0;

mTrajectory = m0;
uTrajectory = [];
wTrajectory = [];

fprintf('the mimial required marking vector: %s\n', mat2str(minMark,4));

K = 0;
stepCnt = 0;
maxStep_last = inf;
maxStep_cur = 0;



% mDif = m - minMark;
% notReach = find(mDif < 0);
% while mDif(notReach)' * mDif(notReach) > 1e-6
while ~isequal((m >= minMark - 1e-6 * ones(PlaceNum,1)), ones(PlaceNum,1))
    % compute the current flows of transitions
    for j = 1:TransNum       
        % enabling degree
        ip = find(Pre(:,j));        
        mi = m(ip);   
        Prei = Pre(:, j);
        PreiNonZero = Prei(ip);
        wmi = mi./PreiNonZero;
        enab(j) = min(wmi);        
        flow(j) = lambda(j) * enab(j);
    end   
        
    stepCnt = stepCnt + 1;
    if stepCnt > K
        % compute the the new fv;
        [w_ss, fv, m_ss, maxStep_cur] = mtof_fv_update(Pre, Post, lambda, delta, m, tj, optFlow);
        %fprintf('current max step: %s, new fv: %s\n', num2str(maxStep_cur, 6), mat2str(fv,4));
        reFv = fv;
        stepCnt = 0;        
    end
    
    if 1
        %****************using modified ON/OFF *****************************%
        for j = 1 : TransNum
            if reFv(j) > 0
                w(j) = min(reFv(j)/delta, flow(j));
                u(j) = flow(j) - w(j);
                reFv(j) = reFv(j) - w(j) * delta;
            else
                ind = find(Pre(:, j));
                if isequal(m(ind) > minMark(ind), ones(size(ind),1) )  
                    difMark = m(ind) - minMark(ind);
                    Preind = Pre(ind, j);
                    w(j) = min(min(difMark ./ (Preind * delta)), flow(j));                
                    u(j) = flow(j) - w(j);   
                else
                    w(j) = 0;
                    u(j) = flow(j);
                end
            end

        end
    else
        %****************using pure ON/OFF%*****************************%
        for j = 1:TransNum 
                w(j) = min(reFv(j)/delta, flow(j));
                u(j) = flow(j) - w(j);
                reFv(j) = reFv(j) - w(j) * delta;
        end
    end
    
    
    
     
    m = m + C * w * delta;
    sumFv = sumFv + w * delta;
    %fprintf('sum fv: %s\n', mat2str(sumFv, 4));
    s = s + 1;
    mTrajectory = [mTrajectory, m];
    wTrajectory = [wTrajectory, w];
    uTrajectory = [uTrajectory, u];
    
%     mDif = m - minMark;
%     notReach = find(mDif < 0);
end

fprintf('steps: %s\n', num2str(s));
%fprintf('required firving count vector: %s\n', mat2str(intialFv, 4));
%fprintf('actural firving count vector: %s\n', mat2str(sumFv, 4));
fprintf('final marking: %s\n', mat2str(m, 4));
