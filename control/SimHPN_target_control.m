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

function SimHPN_target_control(Pre, Post, m0, lambda, sampling, mf, fv)

C = Post - Pre;
[placeNum, transNum] = size(C);

if (size(fv, 1) == 0)
    % Compute minimal firing count vector
    fv = SimHPN_control_minFv(C, m0, mf); 
end

fprintf('fv: %s\n', mat2str(fv,6));

fprintf('choose a control methods: \n 1, Onoff \n 2, onoff+ \n 3, b-onoff \n 4, mpc \n 5, mpc-onoff \n 6, appro.min \n 7, affine \n 8, MPC-plus \n 9, minTime \n');


while(1)
    type = input('chose the control method: ');
    if type == 0
        break;
    elseif type == 1
        [ s, m, u, w ] = SimHPN_control_onoff(Pre, Post, m0, lambda, sampling, mf, fv ); 
        method = 'onoff';
    elseif type == 2
        [ s, m, u, w ] = SimHPN_control_onoff_plus( Pre, Post, m0, lambda, sampling, mf, fv );
        method = 'onoff+';
    elseif type == 3
        d = input('d =');
        [ s, m, u, w ] = SimHPN_control_b_onoff( Pre, Post, m0, lambda, sampling, mf, fv, d ); 
        method = 'b-onoff';
    elseif type == 4
        N = input('N =');
        Q = 1000 * eye(placeNum);
        R = 0.01*eye(transNum);
        %R(1,1) = 0.01;
        Z = SimHPN_control_dare_P( Pre, Post, sampling, Q, R);   
        wf = zeros(size(Pre,2),1);
        [s, m, u, w] = SimHPN_control_mpc( Pre, Post, m0, lambda, sampling, mf, wf, Z, Q, R, N);
        method = 'mpc';
    elseif type == 5
        N = input('N =');
        epsilon = 1e-5;
        zeta = 1e-5;
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
        Z = 1000*eye(placeNum);
        Q = eye(placeNum);
        wf = zeros(size(Pre,2),1);

        [s, m, u, w] = SimHPN_control_mpc_onoff2( Pre, Post, m0, lambda, sampling, mf, fv, wf, Z, Q, R, N, epsilon, zeta);
        method = 'mpc-onoff';
    elseif type == 6
        epsilon = 0.01;
        [time, mTmp, x] = SimHPN_control_approachMin(Pre, Post, lambda, m0, mf, epsilon);
        % get the discrete-time marking in each step
        s = floor(sum(time)/sampling);
        method = 'appro.min';
    elseif type == 7
        [t,m,w,u,F,g] = SimHPN_control_AC_function(Pre,Post,lambda,m0,mf,sampling,0,[1;1]);
        s = size(t, 2);
        method = 'affine';
    elseif type == 8
        N = input('N = ');
        Q = eye(placeNum);
        Z = eye(placeNum) * 1000;
        [s, m, u, w] = SimHPN_control_mpc_plus( Pre, Post, m0, lambda, sampling, mf, fv, N, Z, Q);
        method = 'MPC-plus';
       
    elseif type == 9
        initStep = input('give a intial steps: ');
        [s, m, u, w] = SimHPN_control_minT( Pre, Post, m0, lambda, sampling, mf, initStep );
        method = 'minTime';
    end
    
    fprintf('Method: %s, steps: %d\n', method, s)
    fprintf('choose a control methods: \n 1, Onoff \n 2, onoff+ \n 3, b-onoff \n 4, mpc \n 5, mpc-onoff \n 6, appro.min \n 7, affine \n 8, MPC-plus \n 9, minTime \n');

end
