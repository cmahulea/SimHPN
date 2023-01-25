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
% WARNING: it is only for Marked Graph at this moment 

%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre, Post: Pre- and post- incidence matrix
% m0:        Initial marking
% lambda:    Transitions firing rates
% delta:     Sampling period
% mf:        Final marking
% fv:        The firing count vector, if fv == [], using the minimal
% maxFlow:   The maximal of flow of transitions in optimal steady states 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Output parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% s:  number of steps for reaching the final marking
% m:  |P|*(s+1) matrix, marking trajectory of places, m(:, k) is the
%     state in k_th step, m(:, 1) = m0, m(:, s+1) = mf
% u:  |T|*(s) matrix, control actions
% w:  |T|*(s) matrix, the controlled flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [s, m, u, w] = onoff_maxFlow(Pre, Post, m0, lambda, delta, fv, maxFlow)

% the accumulative amount of firings of transitions
%sumFv = zeros(TransNum, 1);

[PlaceNum, TransNum] = size(Pre);


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
    minMark(i) = maxFlow / lambda(j);
end

m = m0;
s = 0;


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
        
        
        %if j == 3    w(j) = flow(j);     end
    end
    
   
    
    m = m + C * w * delta;
    sumFv = sumFv + w * delta;
    s = s + 1;
end

fprintf('steps: %s\n', num2str(s));
fprintf('minimal required marking: %s\n', mat2str(minMark, 4));
fprintf('required firving count vector: %s\n', mat2str(fv, 4));
fprintf('actural firving count vector: %s\n', mat2str(sumFv, 4));
fprintf('final marking: %s\n', mat2str(m, 4));
