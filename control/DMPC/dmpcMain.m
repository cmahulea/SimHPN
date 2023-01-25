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

function [s, mTrajectory, wTrajectory] = dmpcMain(Pre, Post, lambda, m0, delta, mf, partitions, N, alpha, eps1)

global Pre_G;
global Post_G;

Pre_G = Pre;
Post_G = Post;
                    
Sub_sys = subSysPartition(Pre, Post, m0, mf, lambda, partitions);

numSub = size(Sub_sys, 2);
wGlb = zeros(size(Pre_G,2),1);
m = m0;
s = 0;
eps2 = 1e-6;

pNumSub = zeros(1, numSub);
distSubQd = zeros(1, numSub);

for i = 1 : numSub
    mDif{i} = Sub_sys{i}.m - Sub_sys{i}.mf;
    distSubQd(i) = mDif{i}' * mDif{i};
    
    pNumSub(i) = size(Sub_sys{i}.p,2);   
    Q{i} = eye(pNumSub(i));
    Z{i} = eye(pNumSub(i)) * 1000;
end

subFlag = zeros(1, numSub);
glbFlag = 0;
for i = 1 : numSub
    if distSubQd(i) >= eps2
        subFlag(i) = 1;
    end
    glbFlag = glbFlag + subFlag(i);
end

global mTrajectory;
global wTrajectory;

mTrajectory = m0;
wTrajectory = [];

tStart = tic;

while (glbFlag > 0)
    for i = 1 : numSub
        if subFlag(i) == 0
            Sub_sys{i} = DMPC_step_solver(Sub_sys{i}, delta, 1, Q{i}, Z{i}, 2, alpha, eps1);  
            %fprintf('subsystem %d solve problem 2, at step %d\n', i, s);
        else
            mLast = Sub_sys{i}.m;
            Sub_sys{i} = DMPC_step_solver(Sub_sys{i}, delta, N, Q{i}, Z{i}, 1, alpha, eps1);
            mCurr = Sub_sys{i}.m;
            
            mChange = mCurr - mLast;
            if mChange' * mChange <= eps2   % cannot move to mf
                Sub_sys{i} = DMPC_step_solver(Sub_sys{i}, delta, 1, Q{i}, Z{i}, 2, alpha, eps1);                
            end
        end
        mDif{i} = Sub_sys{i}.m - Sub_sys{i}.mf;
        distSubQd(i) = mDif{i}' * mDif{i};      
        
        %costAll = costAll + FVAL;
        %costAll = costAll + distSubQd(i);
    end
   
    % update flag   
    glbFlag = 0;
    for i = 1 : numSub
        if distSubQd(i) <= eps2
            if subFlag(i) == 1
                fprintf('subsystem %d reached at %d th step\n', i, s+1);            
            end
            subFlag(i) = 0;
        else
            subFlag(i) = 1;
        end
        glbFlag = glbFlag + subFlag(i);
    end
     
      
    %update buffer places states
    for i = 1 : numSub
        for j = 1 : size(Sub_sys{i}.t, 2)
            wGlb(Sub_sys{i}.t(j)) = Sub_sys{i}.w(j);
        end
    end
    m = m + (Post_G - Pre_G) * wGlb * delta;
    wGlb = wGlb;
    

    for i = 1 : numSub
        for biSub = 1 : size(Sub_sys{i}.b, 2)
            biGlb = Sub_sys{i}.b(biSub);
            Sub_sys{i}.bm(biSub) = m(biGlb);
        end
    end
 
    s = s + 1;
    mTrajectory = [mTrajectory  m];
    wTrajectory = [wTrajectory  wGlb];
        
end

telapsed = toc(tStart);

fprintf('id   mf\n')
for i = 1 : size(Pre,1)
    fprintf('%d    ', i);
    if i <= size(Pre,1)
        fprintf('%4.3f', m(i));
    end
    fprintf('\n');
end
fprintf('final states of subsysytems reached in %d steps\n', s);
fprintf('Elapsed time: %f ms\n', telapsed * 1000);

end


