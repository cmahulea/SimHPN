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

% This function is used to compute the makings on the borders of regions in
% the case the flow of transitions are CONSTANT

function [mb, PiRespRowDif] = SimHPN_control_mrkBorder(Pre, Post, m0, mf)

% the colume of mb is the marking on the border
% the row of PiResRowDif is the difference of the rows corresponding to
% the configurations matrix of the both side of the border 

[placeNum, transNum] = size(Pre);
mb = [];
PiRespRowDif = [];

% check every transition to find out the transtions that have multiple
% input places (since these transitions induce the change of
% configuratios). Then build a matrix tPi (colum number equal to placeNum), 
% whose rows to be the rows of cofiguration matrice coresponding to these transitions.

tPiIndex = 1;
mulInputT = 0; % number of transtions with multiple input

for i=1:transNum
    Prei = Pre(:, i);
    nonZero = find(Prei);
    inputNum = size(nonZero, 1);
    if inputNum <= 1
        continue;
    else % t_i has more than one input places
        mulInputT = mulInputT + 1;
        for j = 1:inputNum
            newRow = zeros(1, placeNum);
            newRow(nonZero(j)) = 1/Pre(nonZero(j), i);
            tPi(tPiIndex,:) = newRow;
            tPiIndex = tPiIndex + 1;
        end
        tPiRowIndex(mulInputT) = tPiIndex;     
    end    
end

%m = m0;

mbIndex = 0;

for (i = 1 : mulInputT)
    if i == 1
        startR = 1;
    else
        startR = tPiRowIndex(i - 1);
    end
    endR = tPiRowIndex(i) - 1;

    for (j = startR : (endR -1))
        %fprintf('The transition(%d) with multiple input\n', j)
        for (k = (j + 1) : endR)
            tPiRowA = tPi(j, :);
            tPiRowB = tPi(k, :);
            
            % solve LLP problem to get the points on border
            newMrk = SimHPN_control_mrkBordLpp(tPiRowA-tPiRowB, Pre, Post, m0, mf);
            
            if newMrk == 0
                continue;
            end
            
            % check the two input places with each marking of a transition  
            equalMakPlcIndex = find(tPiRowA-tPiRowB);

            if ((newMrk(equalMakPlcIndex)' * newMrk(equalMakPlcIndex) == 0) || ... 
                    ((newMrk(equalMakPlcIndex) - m0(equalMakPlcIndex))' *  (newMrk(equalMakPlcIndex) - m0(equalMakPlcIndex)) == 0) || ... 
                    ((newMrk(equalMakPlcIndex) - mf(equalMakPlcIndex))' * (newMrk(equalMakPlcIndex) - mf(equalMakPlcIndex)) == 0))
                continue;
            else          
            %if (newMrk' * newMrk ~= 0) & ((newMrk - m0)' *  (newMrk - m0) ~= 0) & ((newMrk - mf)' * (newMrk - mf) ~= 0)
                mbIndex = mbIndex + 1;
                mb(:, mbIndex) = newMrk;
                PiRespRowDif(mbIndex, :) = tPiRowA-tPiRowB; 
            end
            
        end
    end

    
end

% Make the markings in order
[mbR, mbC] = size(mb);

order = 0;
orderMb = 0;
for (i = 1 : placeNum)    
    if (m0(i) > mf(i))
        order = 0; % decrease order
        break;
    elseif m0(i) < mf(i)
        order = 1; % increase order
        break;
    end
end

for i = 1: (mbC -1)
    for j = (i + 1) : mbC
        for (k = 1 : placeNum)
            if mb(k, i) > mb(k, j)
                orderMb = 0; % decrease
                break;
            elseif mb(k, i) < mb(k, j)
                orderMb = 1; % increase
                break;
            end
        end
        
        if orderMb ~= order
            temp = mb(:, i);
            mb(:, i) = mb(:, j);
            mb(:, j) = temp;
        end
    end
end
