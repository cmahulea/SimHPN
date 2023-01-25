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

function [time, m, x] = SimHPN_control_approachMin(Pre, Post, lambda, m0, mf, epsilon)

%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre, Post: Pre- and post- incidence matrix
% m0:        Initial marking
% lambda:    Transitions firing rates
% delta:     Sampling period
% mf:        Final marking
% epsilon:   Stopping threadhold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Output parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time: 1 * (S + 1) matrix, the time costed between two states
%       S is the nomber of board being crossed
% m:  |P|* (S+2) matrix, the system states including m0 and mf, S is the
%     number of intermediate states. m(1) = m0, m(s+2) = mf
% x:  |T|* (S + 1) matrix, accumulative controlled flow (CONSTANT) between states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = [];
x = [];
m = [m0];

global tEscape;
tStart = tic;


% get the marking on the boards that have been crossed
[mb_pw, minTmb_pw, mb_w] = SimHPN_control_mrkBorder_pw(Pre, Post, lambda, m0, mf);
mbNum_pw = size(mb_pw, 2);


if (mbNum_pw == 0)
    %fprintf('no board crossed\n');
    [t, w] = SimHPN_control_minToneR(Pre, Post, lambda, m0, mf);
	[time, interStateNum, x, md] = SimHPN_control_split(Pre, Post, lambda, m0, mf, t, w, epsilon);
    
    m = [m, md, mf];
else
	for (i = 0:mbNum_pw )
		if (i == 0)
			m0_tmp = m0;
			mf_tmp = mb_pw(:, 1);
			t = minTmb_pw(1);
            w = mb_w(:,1);
		elseif (i == mbNum_pw)
			m0_tmp = mb_pw(:, i);
			mf_tmp = mf;
			t =  minTmb_pw(i+1);
            w = mb_w(:, i+1);
		else
			m0_tmp = mb_pw(:, i);
			mf_tmp = mb_pw(:, i+1);
			t = minTmb_pw(i + 1);
            w = mb_w(:, i+1);
        end
    
	[timeTmp, interStateNumtmp, wTmp, mTmp] = SimHPN_control_split(Pre, Post, lambda, m0_tmp, mf_tmp, t, w, epsilon);
    time = [time, timeTmp];
      
    x = [x, wTmp];
    m = [m, mTmp, mf_tmp];
	end
end

% get ending time tick
tEscape = toc(tStart);
fprintf('Time consumed: %f millisecond\n', tEscape * 1000);

 
