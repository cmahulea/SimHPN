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

% SimHPN_control_mpc: standard MPC controller

function [s, m, u, w] = SimHPN_control_mpc( Pre, Post, m0, lambda, delta, mf, wf, Z, Q, R, N)

%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre, Post: Pre- and post- incidence matrix
% m0:        Initial marking
% lambda:    Transitions firing rates
% delta:     Sampling period
% mf:        Final marking
% Z, Q, R:   Matrix
% N:         Horizon step
% fv:        Firing count vector, if not specified, using the minimal one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Output parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% s:  number steps for reaching final state
% m:  |P|* (s+1) matrix, m(:, k) is the final state in k_th step, 
%     m(:,1) = m0, m(:, (S+1)) = mf
% u:  |T|*(s) matrix, control actions
% w:  |T|* s matrix, w(:, k) is the controlled flow in k_th step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off

global tEscape;

tStart = tic;

m = [m0];
u = [];
w = [];

s = 0;
qdis = 1e-6;

temp = 0;

while (1)
    [marking, cflow] = SimHPN_control_mpc_step_QP2( Pre, Post, m0, lambda, delta, mf, wf, Z, Q, R, N );       
    
    s = s + 1;
    
    m_1 = marking(:, 1);
    w_1 = cflow(:, 1);
    
    m = [m, m_1];
    w = [w, w_1];  
    
    dif = m_1 - mf;
    if dif' * dif < qdis  
        fprintf('done!\n');        
        break;
    else
        m0 = m_1;
    end

    temp = temp +1;
    if mod(temp, 100) == 0
        fprintf('-');
    end

end
fprintf('\n');

% get ending time tick
tEscape = toc(tStart);

fprintf('Time consumed: %f millisecond\n', tEscape * 1000);
fprintf('mf can be reached in %d steps (N = %d)\n', s, N);
%Z = Z
%Q = Q
%R = R

warning on



