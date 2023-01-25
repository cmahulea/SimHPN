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

function [is_bound_min,is_bound_max, fs_bound] = SimHPN_Bounds ( pre, post, m0, lambda, control)
% BOUNDS - computes lower and upper bound for infinite servers semnatics
% and upper bound for finite servers semantics AND optimal steady-state for
% controlled contPN system. Only valid for mono-t-semiflow systems.
%
% [is_bound_min,is_bound_max, fs_bound] = bounds ( pre, post, m0, lambda, control)
%       pre     - pre-incidence matrix
%       post    - post-incidence matrix
%       m0      - initial marking
%       lambda  - vector of firing rates
%       control     0: Computes throughput bounds of the steady state.
%                   1: Computes optimal steady state of controlled system.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  RETURNED VALUES FOR THROUGHPUT BOUNDS         %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% is_bound_min : throughput lower bound for transitions
%				  under infinite server semantics
% is_bound_max : throughput upper bound for transitions
%				  under infinite server semantics
% fs_bound : throughput bound for transitions
%				  under finite server semantics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  RETURNED VALUES FOR OPTIMAL STEADY-STATE      %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% is_bound_max(1:nts):                  optimal flow in steady-state
% is_bound_max(nts+1:nts+nps):          marking of optimal flow
% is_bound_max(nts+nps+1:2*nts+nps):    control input for steady-state


nps=size(pre,1);  % Number of places
nts=size(pre,2);  % Number of transitions

% Get p[t]-semiflows
[psemis, tsemis]= SimHPN_ptsemis(post-pre);

if size(tsemis,2)>1
    uiwait(warndlg('Net must be mono-t-semiflow! The following bounds are not valid!'));
    tsemis=tsemis(:,1);
end

% Normalize t-semiflow for transition 1
tsemis=tsemis/tsemis(1);

% IS_BOUND

% This bound is computed by means of the linear programming problem:
%  See "Steady State Performance Evaluation of continuous mono-T-semiflow
%  Petri Nets".
contrtr = [];
contrw = [];
contrz = [];
is_bound_min = [];
if (control == 1)
    temp = mat2str(ones(1,nts));
    temp2 = mat2str(zeros(1,nps));
    contr = char(inputdlg({'Uncontrolled transitions:','Gain vector w.r.t. flow: w''',...
        'Cost vector due to immobilization to maintain the product flow z'''},...
        'Optimal steady-state',1,{'[]',temp,temp2}));
    if isempty(deblank(contr))
        is_bound_max = [];
        is_bound_min = [];
        fs_bound = [];
        return;
    end
    contrtr = setdiff(1:nts,eval(contr(1,:)));
    contrw = eval(contr(2,:));
    contrz = eval(contr(3,:));
    if (size(contrw,1) == 1)
        contrw = contrw';
    end
    if (size(contrz,1) == 1)
        contrz = contrz';
    end
    if (size(contrw,1) ~= nts)
        error('Gain vector w should have dimension equal to the number of transitions!');
    end
    if (size(contrz,1) ~= nps)
        error('Cost vector z should have dimension equal to the number of places!');
    end
end

global bound sol;
eqs = {};
C = post - pre;
ctype = '';
% AA:   f  m sigma
AA=[];
bb=[];

%state equation
AA = [AA;zeros(nps,nts) eye(nps) -C];
bb= [bb;m0];
for i = 1 : nps
    ctype = sprintf('%sS',ctype);
end

for i = 1 : nts %add
    inputPlaces = find(pre(:,i) > 0);
    if (length(inputPlaces) == 1)  %add the equality f(i) = lambda*m/pre
        temp = zeros(1,2*nts+nps);
        temp(i) = 1;
        temp(nts+inputPlaces(1)) = -lambda(i)/pre(inputPlaces(1),i);
        AA = [AA;temp];
        bb = [bb;0];
        ctype = sprintf('%sS',ctype);
    else
        for j = 1 : length(inputPlaces) % for all places add the inequality
            temp = zeros(1,2*nts+nps);
            temp(i) = 1;
            temp(nts+inputPlaces(j)) = -lambda(i)/pre(inputPlaces(j),i);
            AA = [AA;temp];
            bb = [bb;0];
            ctype = sprintf('%sU',ctype);
        end
    end
end

%% equality for a steady state flow
AA = [AA ; C zeros(nps,nps+nts)];
bb = [bb ; zeros(nps,1)];
for i = 1 : nps
    ctype = sprintf('%sS',ctype);
end

vartype = '';
for i = 1 : size(AA,2)
    vartype = sprintf('%sC',vartype);
end


if (control == 0)
    bound = Inf;
    SimHPN_Branch_Bound(1,pre, post, m0, lambda, AA, bb, ctype, vartype, eqs, control, contrtr, contrw, contrz);
    is_bound_min = bound * tsemis;
end

bound = 0;
SimHPN_Branch_Bound(2,pre, post, m0, lambda, AA, bb, ctype, vartype, eqs, control, contrtr, contrw, contrz);
is_bound_max = bound * tsemis;

if (control == 1)
    is_bound_min = sol(1:nts+nps);

    %compute the control input
    mx = sol(nts+1:nts+nps);
    Pre1 = zeros(size(pre,1),size(pre,2));
    for i = 1 : size(pre,2)
        inputP = find(pre(:,i));
        temp = mx(inputP);
        for j = 1 : length(inputP)
            temp(j) = temp(j)/pre(inputP(j),i);
        end
        [mi,ii] = min(temp);
        Pre1(inputP(ii),i) = 1/pre(inputP(ii),i);
    end
    Pre1 = Pre1';
    fxx = diag(lambda)*Pre1*mx;

    is_bound_min = [is_bound_min;fxx-is_bound_max];
    fs_bound = [];
    return;
end


% FS_BOUND

% Bound computed maximizing {k | kv<=lambda}. See "PN Fluidification
%  revisited: Semantics and Steady State". This bound is reached
%  by CF systems.

fs_bound=min(lambda./tsemis)*tsemis;
