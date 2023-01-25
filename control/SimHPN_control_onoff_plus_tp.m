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

% tcpn_onoff() is the implementation of ON-OFF controller for structurally persistent TCPN

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


function [ s, m, u, w ] = SimHPN_control_onoff_plus_tp( Pre, Post, m0, lambda, delta, mf, fv )

% net structure verification 
% ON-OFF controller is only applicable for stucturally persistent TCPN if:
% (A1) The flow matrix C has full rank, OR
% (A2) The net is strongly connected and consistent.

% number of transitions and places

% Norm2 distance threshold to final state
qdis = 1e-6;

tStart = tic;

[placeNum, transNum] = size(Pre);
C = Post - Pre;

% construct matrix G
DT = [];
TA = [];

for (i = 1:placeNum)
    for (j = 1:transNum)
        if (Pre(i, j) > 0)
            tmpD = zeros(1, transNum);
            tmpD(j) = 1;
            tmpTA = zeros(1, placeNum);
            tmpTA(i) = lambda(j)/Pre(i, j);
            DT = [DT;tmpD];
            TA = [TA;tmpTA];
        end
    end
end

G = [DT -TA];

if (size(fv, 1) == 0)
    % Compute minimal firing count vector
    fv = SimHPN_control_minFv(C, m0, mf);    
    fprintf('\nCorresponding minimal firing count vector is: \n');
    fv
else
    fprintf('\nUsing a given firing count vector:\n');
    fv
end

if isequal(fv, zeros(size(Pre,2), 1))
    s = 0;
    return;
end



%%%%%%%%%%%%%%%%%%%Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = [w_{k}, m_{k+1}]'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 0;
m = m0;
u = [];
w = [];
mcur = m0;

MAXPATH = 1;

fv_cur = fv;

fv_cur(1) = 0;
fv_cur(3) = 0;
fv_cur(4) = 0;

f1 = lambda(1) * SimHPN_control_enabling(Pre, 1, mcur);
f2 = lambda(2) * SimHPN_control_enabling(Pre, 2, mcur);
f3 = lambda(3) * SimHPN_control_enabling(Pre, 3, mcur);
f4 = lambda(4) * SimHPN_control_enabling(Pre, 4, mcur);

REF1 = (fv(1) / f1) / (fv(2) / f2);

REF3 = (fv(3) / f3) / (fv(2) / f2);
REF4 = (fv(4) / f4) / (fv(2) / f2);

%count = 0;
%flag = 0;

flag1 = 0;

flag3 = 0;
flag4 = 0;

dif1_cur = REF1
dif3_cur = REF3;
dif4_cur = REF4;


count1 = 0;
count3 = 0;
count4 = 0;

d2 = 10;
fv_cur(1) = fv(1);
flag1 = 1;

f30 = lambda(3) * SimHPN_control_enabling(Pre, 3, m0);
f40 = lambda(3) * SimHPN_control_enabling(Pre, 4, m0);

while ( (mcur - mf)' * (mcur - mf) > qdis )
    mpre = mcur;
    [mcur, ucur, wcur] = SimHPN_control_onoff_plus_LPP( Pre, Post, mpre, lambda, delta, fv_cur);
    m = [m, mcur];
    w = [w, wcur];        
    u = [u, wcur];  % now correct, check later
    fv_cur = fv_cur - wcur * delta;
    
    %d30 = (fv(3)/f3)/(fv(1) / wcur(1))
    %d40 = (fv(4)/f4)/(fv(1) / wcur(1))      
    %halt
%    if flag == 0
%        if count < MAXPATH
            
          
            [mtemp, utemp, wtemp] = SimHPN_control_onoff_plus_LPP( Pre, Post, mpre, lambda, delta, fv_cur); % should use a more specific function to save time
            fv_temp = fv_cur - wtemp * delta;
            
            
 
            
            if flag1 == 0                          
                
                f1 = lambda(1) * SimHPN_control_enabling(Pre, 1, mcur);
                dif1_last = dif1_cur;
                dif1_cur = (fv(1) / f1) / (fv_temp(2) / wtemp(2))
                if (dif1_cur >= dif1_last && count1 >= 3 *2 ) || dif1_cur < d2
                    fv_cur(1) = fv(1);                    
                    s = s
                    flag1 = 1;
                    mcur
                    fv_cur
                else
                    count1 = count1 + 1;
                end
            end
            
            
            if flag3 == 0                       
                
                f3 = lambda(3) * SimHPN_control_enabling(Pre, 3, mcur);
                dif3_last = dif3_cur;
                dif3_cur = (fv(3) / f3) / (fv_temp(2) / wtemp(2))
                if (dif3_cur >= dif3_last && count3 >= 1 ) || dif3_cur < d2
                    fv_cur(3) = fv(3);                    
                    s = s
                    flag3 = 1;
                    mcur
                    fv_cur                    
                else
                    count3 = count3 + 1;
                end
            end
       
            if flag4 == 0 
                f4 = lambda(4) * SimHPN_control_enabling(Pre, 4, mcur);            
                dif4_last = dif4_cur;
                dif4_cur = (fv(4) / f4) / (fv_temp(2) / wtemp(2))
                if (dif4_cur >= dif4_last && count4 >= 1) || dif4_cur < d2
                    s = s
                    fv_cur(4) = fv(4);
                    flag4 = 1;
                    mcur
                    fv_cur
                    s
                else
                    count4 = count4 + 1;
                end
            end

            
            
            %dif_cur = (fv(3) / f3) / (fv(1) / f1)
            %if dif_cur < 3

 %           if dif_cur > dif_last
 %               count = count + 1
 %           end
 %       else
 %           fv_cur(3) = fv(3);
 %           fv_cur(4) = fv(4);
 %           flag = 1;
 %           mcur = mcur
 %           s
 %       end
        
 %   end
    
   
    %{
    if count <= 2 * MAXPATH
        f1 = lambda(1) * SimHPN_control_enabling(Pre, 1, mpre);
        f3 = lambda(3) * SimHPN_control_enabling(Pre, 3, mpre);
        %dif_last = dif_cur;
        dif_cur = (fv(3) / f3) / (fv(1) / f1);
        
        if dif_cur >= REF || dif_cur > dif_last
            count = count + 1
        end
    else
        fv_cur(3) = fv(3);
        fv_cur(4) = fv(4);
    end
    %}

    s = s + 1;
end

step = s
mpre






