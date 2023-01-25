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

% generate subsystem
% Input parameter: Pre, Post, m0, mf, lambda, partitions
%                  partitions{i}.p the places in subsytem i % given
%                  partitions{i}.t the transitions in subsystem i % given
%                  partitions{i}.it the interface transitins in susbystem i
%                  partitions{i}.b the beffer places of subsytem i
% Output parameter: Sub_sys
%                   Sub_sys{i}.Pre the Pre matrix of subsystem i 
%                   Sub_sys{i}.Post the Pre matrix of subsystem i 
%                   Sub_sys{i}.p the places in subsytem i
%                   Sub_sys{i}.t the transitios in subsytem i
%                   Sub_sys{i}.b the buffer places in subsytem i
%                   Sub_sys{i}.bm the state of buffer places
%                   Sub_sys{i}.bmf the final state of buffer places 
%                   Sub_sys{i}.m the state of subsystem i 
%                   Sub_sys{i}.mf the final state of subsystem i 
%                   Sub_sys{i}.lambda the firing rate of subsystem i
%                   Sub_sys{i}.w the flow in each step


function Sub_sys = DCtrl_partion(Pre, Post, m0, mf, lambda, partitions)

[row, column] = size(partitions);
num_sub = column; 

for i = 1 : num_sub
    pIdx = partitions{i}.p;
    tIdx = partitions{i}.t;
    bInf = partitions{i}.b;
    
    for pn = 1 : size(pIdx,2)
        Sub_sys{i}.m(pn,1) = m0(pIdx(pn),1); 
        Sub_sys{i}.mf(pn,1) = mf(pIdx(pn),1); 
        for tn = 1 : size(tIdx,2)
            Sub_sys{i}.Pre(pn,tn) = Pre(pIdx(pn), tIdx(tn));
            Sub_sys{i}.Post(pn,tn) = Post(pIdx(pn), tIdx(tn));
            Sub_sys{i}.lambda(tn) = lambda(tIdx(tn));
        end
    end
    
    for bn = 1 : size(bInf,2)
        Sub_sys{i}.bm(bn) = m0(bInf(1, bn));
        Sub_sys{i}.bmf(bn) = mf(bInf(1, bn));
    end
    
    Sub_sys{i}.p = pIdx;    
    Sub_sys{i}.t = tIdx;    
    Sub_sys{i}.b = bInf;     
    Sub_sys{i}.it = partitions{i}.it; 
    Sub_sys{i}.w = [];
end

% for i = 1 : num_sub
%     fprintf('subsystem %s:\n', num2str(i));
%     fprintf('.p: %s\n', mat2str(Sub_sys{i}.p));
%     fprintf('.t: %s\n', mat2str(Sub_sys{i}.t));
%     fprintf('.lambda: %s\n', mat2str(Sub_sys{i}.lambda));
%     fprintf('.bInf: %s\n', mat2str(Sub_sys{i}.b));
%     fprintf('.Pre: %s\n', mat2str(Sub_sys{i}.Pre));
%     fprintf('.Post: %s\n', mat2str(Sub_sys{i}.Post));
%     fprintf('.m: %s\n', mat2str(Sub_sys{i}.m));
%     fprintf('.mf: %s\n', mat2str(Sub_sys{i}.mf));
%     fprintf('.bm: %s\n', mat2str(Sub_sys{i}.bm));
%     fprintf('.bmf: %s\n', mat2str(Sub_sys{i}.bmf));
%     fprintf('.w: %s\n', mat2str(Sub_sys{i}.w));
%     
%     fprintf('-------------------------------------\n\n');
% end


