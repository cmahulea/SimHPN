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

% SimHPN_control_law(): main function, compute the control law to reach a
% given final state mf from m0

%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre, Post: Pre- and post- incidence matrix
% lambda:    Transitions firing rates
% sampling:  Sampling period
% m0:        Initial marking
% method:       control methods
%            'onoff': standard ON/OFF controller (only for CF net system)
%            'onoff_plus': propotional ON/OFF controller
%            'b_onoff_plus': balanced ON/OFF controller
%            'mpc': standard MPC controller
%            'mpc_onoff': MPC-ON/OFF controller
%            'appMin': approaching mimum-time controller (Hanife's paper)
%            'affCtrl': affine controller
%            'minT': mininum-time controller (may be very expensive)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Output parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
% s:  number of steps for reaching the final marking
% m:  |P|*(s+1) matrix, marking trajectory of places, m(:, k) is the
%     state in k_th step, m(:, 1) = m0, m(:, s+1) = mf
% u:  |T|*(s) matrix, control actions
% w:  |T|*(s) matrix, the controlled flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ s, m, u, w ] = SimHPN_control_law( Pre , Post , lambda , m0 , sampling, method )

[placeNum, transNum] = size(Pre);
        
mftemp = []; % fix the final state, for test use.
        
switch method
    case 'onoff'
        prompt={'The final state Mf:','Firing count vector (if not specified, using the minimal one):'};
        name='ON-OFF Controller';
        numlines=1;
        default={mat2str(mftemp),'[]'};
        options.Resize='on';
        options.WindowStyle='normal';
        
        contr = char(inputdlg(prompt,name,numlines,default, options));
        
        if isempty(deblank(contr))
            mf = [];
            fv = [];
        else
            mf = eval(contr(1,:));
            fv = eval(contr(2,:));        
        end
                
        inputChk = 1;
        if isempty(mf)
            inputChk = 0;
        elseif ~isequal(size(mf), size(m0)) 
            errordlg('The dimension of Mf is not correct','SimHPN: Data error');
            inputChk = 0;
        elseif ~isempty(fv) && ~isequal(size(fv), size(lambda)) 
            errordlg('The dimension of firing count vector is not correct','SimHPN: Data error');
            inputChk = 0;
        end
  
        if inputChk == 0
            s = 0;
            m = [];
            u = [];
            w = [];
            return;
        end

        [ s, m, u, w ] = SimHPN_control_onoff(Pre, Post, m0, lambda, sampling, mf, fv );     
       
        if s == 0
            m = [];
            u = [];
            w = [];
            return;            
        end        
    case 'onoff_plus'
        prompt={'The final state Mf:','Firing count vector (if not specified, using a minimal one):'};
        name='ON-OFF-plus Controller';
        numlines=1;
        default={mat2str(mftemp),'[]'};
        options.Resize='on';
        options.WindowStyle='normal';
        
        contr = char(inputdlg(prompt,name,numlines,default, options));
        
        if isempty(deblank(contr))
            mf = [];
            fv = [];
        else
            mf = eval(contr(1,:));
            fv = eval(contr(2,:));        
        end
                
        inputChk = 1;
        if isempty(mf)
            inputChk = 0;
        elseif ~isequal(size(mf), size(m0)) 
            errordlg('The dimension of Mf is not correct','SimHPN: Data error');
            inputChk = 0;
        elseif ~isempty(fv) && ~isequal(size(fv), size(lambda)) 
            errordlg('The dimension of firing count vector is not correct','SimHPN: Data error');
            inputChk = 0;
        end
  
        if inputChk == 0
            s = 0;
            m = [];
            u = [];
            w = [];
            return;
        end
        [ s, m, u, w ] = SimHPN_control_onoff_plus( Pre, Post, m0, lambda, sampling, mf, fv );
       
        if s == 0
            m = [];
            u = [];
            w = [];
            return;            
        end

        if size(u,1) == 0
            % computer the control action u
            for step = 1 : s
                curM = m(:, step);
                for j = 1 : transNum
                    ucFlow(j, step) = SimHPN_control_enabling(Pre, j, curM) * lambda(j);
                end
            end
            %t = ones(1, s) * delta;
            u = ucFlow -  w;
        end        
        
    case 'b_onoff'
        prompt={'The final state Mf:','Firing count vector (if not specified, using a minimal one):', 'd:'};
        name='Balanced ON-OFF-plus Controller';
        numlines=1;
        default={mat2str(mftemp),'[]', '10'};
        options.Resize='on';
        options.WindowStyle='normal';
        
        contr = char(inputdlg(prompt,name,numlines,default, options));
        
        if isempty(deblank(contr))
            mf = [];
            fv = [];
        else
            mf = eval(contr(1,:));
            fv = eval(contr(2,:));
            d = eval(contr(3,:));
        end
                
        inputChk = 1;
        if isempty(mf)
            inputChk = 0;
        elseif ~isequal(size(mf), size(m0)) 
            errordlg('The dimension of Mf is not correct','SimHPN: Data error');
            inputChk = 0;
        elseif ~isempty(fv) && ~isequal(size(fv), size(lambda)) 
            errordlg('The dimension of firing count vector is not correct','SimHPN: Data error');
            inputChk = 0;
        elseif d < 1
            errordlg('d should be greater than 1','SimHPN: Data error');
            inputChk = 0;
        end
  
        if inputChk == 0
            s = 0;
            m = [];
            u = [];
            w = [];
            return;
        end

        [ s, m, u, w ] = SimHPN_control_b_onoff( Pre, Post, m0, lambda, sampling, mf, fv, d );        
       
        if s == 0
            m = [];
            u = [];
            w = [];
            return;            
        end

        if size(u,1) == 0
            % computer the control action u
            for step = 1 : s
                curM = m(:, step);
                for j = 1 : transNum
                    ucFlow(j, step) = SimHPN_control_enabling(Pre, j, curM) * lambda(j);
                end
            end
            %t = ones(1, s) * delta;
            u = ucFlow -  w;
        end          
        
    case 'mpc'
        prompt={'The final state Mf:', 'The final contrlled flow: Wf', 'Horizon step N:', ... 
            'Matrix Z:', 'Matrix Q:', 'Matrix R:'};
        name='MPC Controller';
        numlines=1;
        

        wf = zeros(size(Pre,2),1);  % fix the final flow, for test use
        
        default={mat2str(mftemp), mat2str(wf), '1', mat2str(1000 * eye(placeNum)), mat2str(1000 * eye(placeNum)), mat2str(0.01 * eye(transNum))};
        options.Resize='on';
        options.WindowStyle='normal';
        
        contr = char(inputdlg(prompt,name,numlines,default, options));

        if isempty(deblank(contr))
            mf = [];
        else
            mf = eval(contr(1,:));
            wf = eval(contr(2,:));
            N = eval(contr(3,:));
            Z = eval(contr(4,:));
            Q = eval(contr(5,:));
            R = eval(contr(6,:));
        end

        inputChk = 1;       
        if isempty(mf)
            inputChk = 0;
        elseif ~isequal(size(mf), size(m0))
            errordlg('The dimension of Mf is not correct','SimHPN: Data error');
            inputChk = 0;
        elseif N < 1 || ~(rem(N,1) == 0)
            errordlg('Horizon step should be a nature number','SimHPN: Data error');
            inputChk = 0;
        elseif ~isequal(size(Z), size(eye(placeNum))) || ~isequal(size(Q), size(eye(placeNum))) || ~isequal(size(R), size(eye(transNum)))
                errordlg('One of the dimensions of Z, Q, R is not correct','SimHPN: Data error');
                inputChk = 0;            
        end       
  
        if inputChk == 0
            s = 0;
            m = [];
            u = [];
            w = [];
            return;
        end
                
        Z = SimHPN_control_dare_P( Pre, Post, sampling, Q, R);   
             
        [s, m, u, w] = SimHPN_control_mpc( Pre, Post, m0, lambda, sampling, mf, wf, Z, Q, R, N);
        
        if s == 0
            m = [];
            u = [];
            w = [];
            return;            
        end
        
        if size(u,1) == 0
            % computer the control action u
            for step = 1 : s
                curM = m(:, step);
                for j = 1 : transNum
                    ucFlow(j, step) = SimHPN_control_enabling(Pre, j, curM) * lambda(j);
                end
            end
            %t = ones(1, s) * delta;
            u = ucFlow -  w;
        end
        
%     case 'mpc_onoff'
%         prompt={'The final state Mf:', 'The final contrlled flow: Wf', 'Horizon step N:', ...
%             'Firing count vector (if not specified, using one of the minimal):', 'Matrix Z:', 'Matrix Q:', 'Matrix R:'};
%         name='MPC Controller';
%         numlines=1;
%         
%         wf = zeros(size(Pre,2),1);
%         
%         default={mat2str(mftemp), mat2str(wf), '1', '[]', mat2str(1000 * eye(placeNum)), ... 
%             mat2str(1000 * eye(placeNum)), mat2str(0.01 * eye(transNum))};
%         options.Resize='on';
%         options.WindowStyle='normal';
%         
%         contr = char(inputdlg(prompt,name,numlines,default, options));
% 
%         if isempty(deblank(contr))
%             mf = [];
%         else
%             mf = eval(contr(1,:));
%             wf = eval(contr(2,:));
%             N = eval(contr(3,:));
%             fv = eval(contr(4,:));            
%             Z = eval(contr(5,:));
%             Q = eval(contr(6,:));
%             R = eval(contr(7,:));
%         end
% 
%         inputChk = 1;       
%         if isempty(mf)
%             inputChk = 0;
%         elseif ~isequal(size(mf), size(m0))
%             errordlg('The dimension of Mf is not correct','SimHPN: Data error');
%             inputChk = 0;
%         elseif N < 1 || ~(rem(N,1) == 0)
%             errordlg('Horizon step should be a nature number','SimHPN: Data error');
%             inputChk = 0;
%         elseif ~isequal(size(Z), size(eye(placeNum))) || ~isequal(size(Q), size(eye(placeNum))) || ~isequal(size(R), size(eye(transNum)))
%                 errordlg('One of the dimensions of Z, Q, R is not correct','SimHPN: Data error');
%                 inputChk = 0;            
%         end       
%   
%         if inputChk == 0
%             s = 0;
%             m = [];
%             u = [];
%             w = [];
%             return;
%         end
%         
%         Z = SimHPN_control_dare_P( Pre, Post, sampling, Q, R);
%              
%         [s, m, u, w] = SimHPN_control_mpc_onoff( Pre, Post, m0, lambda, sampling, mf, fv, wf, Z, Q, R, N);
%         
%         if s == 0
%             m = [];
%             u = [];
%             w = [];
%             return;            
%         end
%         
%         if size(u,1) == 0
%             % computer the control action u
%             for step = 1 : s
%                 curM = m(:, step);
%                 for j = 1 : transNum
%                     ucFlow(j, step) = SimHPN_control_enabling(Pre, j, curM) * lambda(j);
%                 end
%             end
%             %t = ones(1, s) * delta;
%             u = ucFlow -  w;        
%         end
    case 'mpc_onoff'
        prompt={'The final state Mf:', 'The final contrlled flow: Wf', 'Time horizon N:', ... 
            'Firing count vector (if not specified, using a minimal one):', ...
            'Matrix Q:', 'Matrix R:'};
        name='MPC-ON/OFF Controller';
        numlines=1;
        
                        
        % compute the matrix R        
        conflictT = zeros(transNum, transNum); % temp(j,j) = 1, if tj is in a conflit
        for i = 1 : placeNum
            row = Pre(i,:);
            if sum(row>0) > 1
                for j = 1 : transNum
                    if Pre(i,j) > 0
                        conflictT(j, j) = 1;
                    end
                end
            end
        end
        
        % w_k(tj) = min{sigma_remain(tj), flow(tj)}, for persistent tj --- ON-OFF
        R = zeros(transNum, transNum);
        for j = 1 : transNum
            R(j, j) = 1000 * (conflictT(j,j) == 0);
        end
                              
        default={mat2str(mftemp), mat2str(zeros(transNum,1)), '1', '[]', mat2str(eye(placeNum)), mat2str(R)};
        options.Resize='on';
        options.WindowStyle='normal';
        
        contr = char(inputdlg(prompt,name,numlines,default, options));

        if isempty(deblank(contr))
            mf = [];
        else
            mf = eval(contr(1,:));
            wf = eval(contr(2,:));
            N = eval(contr(3,:));
            fv = eval(contr(4,:));            
            %Z = eval(contr(5,:));
            Q = eval(contr(5,:));
            Z = Q;
            R = eval(contr(6,:));
            %epsilon = eval(contr(8,:));
            %zeta = eval(contr(9,:));
        end
        
        epsilon = 1e-5;
        zeta = 1e-5;
                

        inputChk = 1;       
        if isempty(mf)
            inputChk = 0;
        elseif ~isequal(size(mf), size(m0))
            errordlg('The dimension of Mf is not correct','SimHPN: Data error');
            inputChk = 0;
        elseif N < 1 || ~(rem(N,1) == 0)
            errordlg('Time horizon should be a nature number','SimHPN: Data error');
            inputChk = 0;
        elseif ~isequal(size(Z), size(eye(placeNum))) || ~isequal(size(Q), size(eye(placeNum))) || ~isequal(size(R), size(eye(transNum)))
                errordlg('One of the dimensions of Z, Q, R is not correct','SimHPN: Data error');
                inputChk = 0;            
        end       
  
        if inputChk == 0
            s = 0;
            m = [];
            u = [];
            w = [];
            return;
        end
        
           [s, m, u, w] = SimHPN_control_mpc_onoff2( Pre, Post, m0, lambda, sampling, mf, fv, wf, Z, Q, R, N, epsilon, zeta);
        
        if s == 0
            m = [];
            u = [];
            w = [];
            return;            
        end
        
        if size(u,1) == 0
            % computer the control action u
            for step = 1 : s
                curM = m(:, step);
                for j = 1 : transNum
                    ucFlow(j, step) = SimHPN_control_enabling(Pre, j, curM) * lambda(j);
                end
            end
            %t = ones(1, s) * delta;
            u = ucFlow -  w;  
        end
    case 'appMin'
        prompt={'The final state Mf:','Epsilon (stopping threshold):'};
        name='Approaching Mininum-time Controller';
        numlines=1;
        default={'[]','0.01'};
        options.Resize='on';
        options.WindowStyle='normal';
        
        contr = char(inputdlg(prompt,name,numlines,default, options));
        
        if isempty(deblank(contr))
            mf = [];
        else
            mf = eval(contr(1,:));
            epsilon = eval(contr(2,:));       
        end
        
        inputChk = 1;       
        if isempty(mf)
            inputChk = 0;       
        elseif ~isequal(size(mf), size(m0))
            errordlg('The dimension of Mf is not correct','SimHPN: Data error');
            inputChk = 0;
        elseif ~isempty(find(m0==0, 1)) || ~isempty(find(mf==0, 1))
            errordlg('This control method can only be applied for possitive initial and final marking','SimHPN: Data error');
            inputChk = 0;
        elseif epsilon >= 1 || epsilon <= 0
            errordlg('Please choose 0 < epsilon < 1','SimHPN: Data error');
            inputChk = 0;
        end
        
        if inputChk == 0
            s = 0;
            m = [];
            u = [];
            w = [];
            return;
        end
        
        [time, mTmp, x] = SimHPN_control_approachMin(Pre, Post, lambda, m0, mf, epsilon);
        
        % get the discrete-time marking in each step
        s = floor(sum(time)/sampling);
        
        sumT = 0;
        timeAcumu = zeros(1, size(time,2));
        for idx = 1 : size(time,2)
            wTmp(:, idx) = x(:, idx) / time(idx);     % controlled flow
            timeAcumu(1,idx) = sumT + time(1,idx);
            sumT = timeAcumu(1,idx);
        end
                        
        m = [m0];        
        for step = 1 : s
            tau = step * sampling;
            for idx = 1 : size(timeAcumu, 2)
                if tau <= timeAcumu(idx)
                    w(:, step) = wTmp(:, idx);
                    break;
                end
            end
            if idx == 1
                lastT = 0;
            else
                lastT = timeAcumu(idx - 1);
            end
            lastM = mTmp(:, idx);
            nextM = mTmp(:, idx+1);
            nextT = timeAcumu(idx);
                           
            curM = lastM + ((nextM - lastM) / time(idx)) * (tau - lastT);                       
            m = [m, curM];
        end
        
               % computer the control action u
        for step = 1 : s
            curM = m(:, step);
            for j = 1 : transNum
                ucFlow(j, step) = SimHPN_control_enabling(Pre, j, curM) * lambda(j);
            end
        end
        %t = ones(1, s) * delta;
        u = ucFlow -  w;                    
    
    case 'affCtrl'
        prompt={'The final state Mf:'};
        name='Affine Controller';
        numlines=1;
        default={'[]'};
        options.Resize='on';
        options.WindowStyle='normal';
        
        contr = char(inputdlg(prompt,name,numlines,default, options));
        
        if isempty(deblank(contr))
            mf = [];
        else
            mf = eval(contr(1,:));
        end
        
        inputChk = 1;       
        if isempty(mf)
            inputChk = 0;       
        elseif ~isequal(size(mf), size(m0))
            errordlg('The dimension of Mf is not correct','SimHPN: Data error');
            inputChk = 0;
        elseif ~isempty(find(m0==0, 1)) || ~isempty(find(mf==0, 1))
            errordlg('This control method can only be applied for possitive initial and final marking','SimHPN: Data error');
            inputChk = 0;
        end
        
        if inputChk == 0
            s = 0;
            m = [];
            u = [];
            w = [];
            return;
        end
        
        [t,m,w,u,F,g] = SimHPN_control_AC_function(Pre,Post,lambda,m0,mf,sampling,0,[1;1]);
        s = size(t, 2);
    case 'minT'
        prompt={'The final state Mf:', 'Initial guess of number of time steps:'};
        name='The mimimum-time control laws';
        numlines=1;
        default={'[]', '0'};
        options.Resize='on';
        options.WindowStyle='normal';
        
        contr = char(inputdlg(prompt,name,numlines,default, options));
        
        if isempty(deblank(contr))
            mf = [];
        else
            mf = eval(contr(1,:));
            initStep = eval(contr(2,:));
        end
        
        inputChk = 1;       
        if isempty(mf)
            inputChk = 0;       
        elseif ~isequal(size(initStep), [1,1]) || initStep <= 0 || mod(initStep,1) ~= 0
            errordlg('Initial steps number should be a nature number','SimHPN: Data error');
            inputChk = 0;
        end
        
         if inputChk == 0
            s = 0;
            m = [];
            u = [];
            w = [];
            return;
        end       
        
        [s, m, u, w] = SimHPN_control_minT( Pre, Post, m0, lambda, sampling, mf, initStep );
        
        if s == 0
            s = 0;
            m = [];
            u = [];
            w = [];
            return;
        end
        
        if size(u,1) == 0
            % compute the control action u
            for step = 1 : s
                curM = m(:, step);
                for j = 1 : transNum
                    ucFlow(j, step) = SimHPN_control_enabling(Pre, j, curM) * lambda(j);
                end
            end
            %t = ones(1, s) * delta;
            u = ucFlow -  w; 
        end
    case 'dmpc'
        prompt={'Number of subsystems:'};
        name='Distributed MPC control';
        numlines=1;
        default={'2'};
        options.Resize='on';
        options.WindowStyle='normal';
        
        contr = char(inputdlg(prompt,name,numlines,default, options));
        
        numSub = 0;
        
        if isempty(deblank(contr))
            s = 0;
            u = [];
            w = [];
            m = [];
            return
        else
            numSub = eval(contr(1,:));
        end
        
        if numSub > 1
            
            for i = 1 : numSub
                temp = strcat('Subsystem ', num2str(i), ':');
                prompt{i}= temp;
                default{i} = ['[]'; '[]'];
            end
            
            prompt{i+1}= 'The final state Mf:';
            default{i+1} = '[]';
            prompt{i+2}= 'Time Horizon N:';
            default{i+2} = '1';
            prompt{i+3}= 'alpha:';
            default{i+3} = '0.5';
            prompt{i+4}= 'epsilon:';
            default{i+4} = '0.1';          
            
            name='Distributed MPC control';
            numlines=2;          
            options.Resize='on';
            options.WindowStyle='normal';
            
            contr = char(inputdlg(prompt,name,numlines,default, options));      
            
            if isempty(deblank(contr))
                s = 0;
                u = [];
                w = [];
                m = [];
                return;
            end 
            
             for i = 1 : numSub
                  partitions{i}.b = eval(contr(2*i-1,:));
                  partitions{i}.it = eval(contr(2*i,:));

             end                      
                  
            mf = eval(contr(2*i+1,:));
            N = eval(contr(2*i+2,:));
            alpha = eval(contr(2*i+3,:));
            eps1 = eval(contr(2*i+4,:));
            
           

            [s, m, w] = dmpcMain(Pre, Post, lambda, m0, sampling, mf, partitions, N, alpha, eps1);

            % computer the control action u
            for step = 1 : s
                curM = m(:, step);
                for j = 1 : transNum
                    ucFlow(j, step) = SimHPN_control_enabling(Pre, j, curM) * lambda(j);
                end
            end
            %t = ones(1, s) * delta;
            u = ucFlow -  w;  
        
        else 
            errordlg('At least 2 subsystems!','SimHPN: Data error');
            s = 0;
            u = [];
            w = [];
            m = [];
        end
   case 'mfc'
       tj = 1;
       tIdxToMax = 0 * ones(1, size(Pre,2));
       tIdxToMax(tj) = 1;
       [sigmaV, m_steady, flow_steady] = optCtl(Pre, Post, lambda, m0, tIdxToMax, zeros(1, size(Pre,1)));
       flowMax = flow_steady(tj);
       [mfv, m, flow] = minFv_optCtl(Pre, Post, lambda, m0, tj, flowMax);
       [s, m, u, w] = onoff_maxFlow_cf(Pre, Post, m0, lambda, sampling, mfv, tj, flowMax);
       fprintf('maximal flow in steady state: %s\n', mat2str(flow_steady, 4));
       
end
