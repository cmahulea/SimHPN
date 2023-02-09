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


function ret = SimHPN(action,varargin)

% Simulator for hybrid petri net systems

% actions:

%  ini:  figure and controls
%  load: load petri net from a .rdp file
%  simula: simulate petri net system
%  ptsemis: compute p t - semiflows of the net
%  bounds: compute throughput bounds (finite and infinite servers
%  semantics)
%             only valid for mono-t-semiflow

if nargin<1,
    action='ini';
end;


thisfile='SimHPN';

switch action

    %================
    % INITIALIZATION
    %================
    case 'ini'  % Initialize figure and controls
        
%         addpath(sprintf('%s/%s',temp,'glpkmex'));
%         addpath(sprintf('%s/%s',temp,'control'));
%         addpath(sprintf('%s/%s',temp,'cdd'));
        
        figpri=figure( ...
            'Name','SimHPN - Hybrid Petri Nets Simulator', ...
            'Position',[35 50 900 600], ...
            'NumberTitle','off', ...
            'MenuBar','none',...
            'ToolBar','auto',...
            'Visible','off'...
            );
        ret = figpri;
        fromx=0.05;tox=0.79;
        fromy=0.4;toy=0.95;

        axes('position',[fromx fromy tox-fromx toy-fromy],'Color',[0 0 0]);
        a = uimenu('Label','Model');
        uimenu(a,'Label','Import from &Pmeditor','Callback',strcat(thisfile,'(''load_rdp'')'));
        uimenu(a,'Label','Import from &TimeNet','Callback',strcat(thisfile,'(''load_tn'')'));
        uimenu(a,'Label','Import from .&mat file','Callback',strcat(thisfile,'(''load_mat'')'));
        uimenu(a,'Label','Save to .&mat file','Callback',strcat(thisfile,'(''save_mat'')'));

        a = uimenu('Label','Options');
        uimenu(a,'Label','Show Figure Toolbar','Callback',strcat(thisfile,'(''opttool'')'), ...
            'Checked','off');
        

        a = uimenu('Label','Simulation');
        uimenu(a,'Label','&Markings to plot','Callback',strcat(thisfile,'(''plot_var_m'')'));
        uimenu(a,'Label','&Flows to plot','Callback',strcat(thisfile,'(''plot_var_f'')'));
        uimenu(a,'Label','&Save results to workspace','Callback',strcat(thisfile,'(''save_simul'')'));

        a = uimenu('Label','Structural');
        uimenu(a,'Label','Conservativeness','Callback',strcat(thisfile,'(''conservativeness'')'));
        uimenu(a,'Label','Consistency','Callback',strcat(thisfile,'(''consistency'')'));
        uimenu(a,'Label','P T semiflows','Callback',strcat(thisfile,'(''ptsemis'')'));

                
                
        a = uimenu('Label','Discrete');
        
        a = uimenu('Label','Continuous');
        b = uimenu(a,'Label','Compute bounds','Callback',strcat(thisfile,'(''bounds'')'));

        b = uimenu(a,'Label','Centralized control');   
        
        menuCtrOnoff = uimenu(b,'Label','&ON/OFF Methods');
        uimenu(menuCtrOnoff,'Label','Standard ON/OFF','Callback',strcat(thisfile,'(''contr_law'',''onoff'')'));
        uimenu(menuCtrOnoff,'Label','ON/OFF-plus','Callback',strcat(thisfile,'(''contr_law'',''onoff_plus'')'));
        uimenu(menuCtrOnoff,'Label','Balanced ON/OFF','Callback',strcat(thisfile,'(''contr_law'',''b_onoff'')'));
        uimenu(menuCtrOnoff,'Label','MPC-ON/OFF','Callback',strcat(thisfile,'(''contr_law'',''mpc_onoff'')'));
        %uimenu(menuCtrOnoff,'Label','MPC-ON-OFF2','Callback',strcat(thisfile,'(''contr_law'',''mpc_onoff2'')'));
        
        uimenu(b,'Label','MPC','Callback',strcat(thisfile,'(''contr_law'',''mpc'')'));
        uimenu(b,'Label','Approaching Min. Time','Callback',strcat(thisfile,'(''contr_law'',''appMin'')'));
        uimenu(b,'Label','Affine Control','Callback',strcat(thisfile,'(''contr_law'',''affCtrl'')'));
        uimenu(b,'Label','Mininum-time Control Laws','Callback',strcat(thisfile,'(''contr_law'',''minT'')'));        
        
        
               
%        uimenu(b,'Label','&Control law','Callback',strcat(thisfile,'(''contr_law'')'));
%        uimenu(b,'Label','&Save results to workspace','Callback',strcat(thisfile,'(''save_control'')'));

        b = uimenu(a,'Label','Distributed control', 'Callback',strcat(thisfile,'(''contr_law'',''dmpc'')')); 
        %uimenu(b,'Label','Decentralized minimum-time for choice-free nets','Callback',strcat(thisfile,'(''contr_law'',''dcenCF'')'));
        %uimenu(b,'Label','Distributed MPC','Callback',strcat(thisfile,'(''contr_law'',''dmpc'')'));
        
        b = uimenu(a,'Label','Minimum-time flow control (CF)', 'Callback',strcat(thisfile,'(''contr_law'',''mfc'')')); 

        b = uimenu(a,'Label','&Save results to workspace','Callback',strcat(thisfile,'(''save_control'')')); 
         
        b = uimenu(a,'Label','Optimal');
        uimenu(b,'Label','Optimal &Observability','Callback',strcat(thisfile,'(''optobs'')'));
        uimenu(b,'Label','Optimal &Control','Callback',strcat(thisfile,'(''optsteady'')'));
        
        b = uimenu(a,'Label','Structural controllability analysis');
        uimenu(b,'Label','Influence of controllable transitions','Callback',strcat(thisfile,'(''influ'')'));
        uimenu(b,'Label','Net &rank-controllability test','Callback',strcat(thisfile,'(''nrc'')'));
        
        b = uimenu(a,'Label','Diagnosis','Callback',strcat(thisfile,'(''diagnosis'')'));

        
        
        a = uimenu('Label','Hybrid');


        % Right frame
        fromx=0.8; tox=0.98;
        fromy=0.02;toy=0.98;

        uicontrol( ...
            'Style','frame', ...
            'Units','normalized', ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'BackgroundColor',[0.70 0.70 0.70]);

        % Controls on the right frame

        % Firing Semantics
        fromx=0.81; tox=0.97;
        fromy=0.8;toy=0.97;
        uicontrol( ...
            'Style','frame', ...
            'Units','normalized', ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'BackgroundColor',[0.80 0.80 0.80]);

        fromx=0.82; tox=0.96;
        fromy=0.91;toy=0.96;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','Firing Semantics:');

        fromx=0.82; tox=0.96;
        fromy=0.86;toy=0.91;
        uicontrol( ...
            'Style','radiobutton', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Value', 1, ...
            'Tag', 'infser', ...
            'CallBack',strcat(thisfile,'(''infser'')'), ...
            'String','Infinite Servers');

        fromx=0.82; tox=0.96;
        fromy=0.81;toy=0.86;

        uicontrol( ...
            'Style','radiobutton', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Tag', 'prod', ...
            'CallBack',strcat(thisfile,'(''prod'')'), ...
            'String','Product Enabling');


        % Plot
        fromx=0.81; tox=0.97;
        fromy=0.57;toy=0.79;
        uicontrol( ...
            'Style','frame', ...
            'Units','normalized', ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'BackgroundColor',[0.80 0.80 0.80]);

        fromx=0.82; tox=0.96;
        fromy=0.73;toy=0.78;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','Plot:');

        fromx=0.82; tox=0.96;
        fromy=0.68;toy=0.73;
        uicontrol( ...
            'Style','radiobutton', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Tag','pmarking', ...
            'Value', 1, ...
            'CallBack',strcat(thisfile,'(''plot_marking'')'), ...
            'String','Marking Evolution');

        fromx=0.82; tox=0.96;
        fromy=0.63;toy=0.68;
        uicontrol( ...
            'Style','radiobutton', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Tag','pflow', ...
            'Value', 0, ...
            'CallBack',strcat(thisfile,'(''plot_flow'')'), ...
            'String','Flow Evolution');

        fromx=0.82; tox=0.84;
        fromy=0.58;toy=0.63;
        uicontrol( ...
            'Style','radiobutton', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Tag','pplaces', ...
            'Value', 0, ...
            'CallBack',strcat(thisfile,'(''plot_mxm'')'), ...
            'Position',[fromx fromy tox-fromx toy-fromy]);

        fromx=0.84; tox=0.86;
        fromy=0.58;toy=0.63;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','m(');

        fromx=0.86; tox=0.89;
        fromy=0.58;toy=0.63;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Tag', 'px', ...
            'String','1');

        fromx=0.89; tox=0.92;
        fromy=0.58;toy=0.63;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','),m(');

        fromx=0.92; tox=0.95;
        fromy=0.58;toy=0.63;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Tag', 'py', ...
            'String','2');

        fromx=0.95; tox=0.96;
        fromy=0.58;toy=0.63;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String',')');

        % Errors

        fromx=0.81; tox=0.97;
        fromy=0.34;toy=0.56;
        uicontrol( ...
            'Style','frame', ...
            'Units','normalized', ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'BackgroundColor',[0.80 0.80 0.80]);

        fromx=0.82; tox=0.96;
        fromy=0.52;toy=0.55;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','Absolute Error:');

        fromx=0.83; tox=0.95;
        fromy=0.48;toy=0.52;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Tag', 'eabs', ...
            'String','1e-6');

        fromx=0.82; tox=0.96;
        fromy=0.45;toy=0.48;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','Relative Error:');

        fromx=0.83; tox=0.95;
        fromy=0.41;toy=0.45;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Tag', 'erel', ...
            'String','1e-3');

        fromx=0.82; tox=0.96;
        fromy=0.38;toy=0.41;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','Sampling Time:');

        fromx=0.83; tox=0.95;
        fromy=0.345;toy=0.385;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Tag', 'sampling', ...
            'String','0.01');

        % Simulate button

%         fromx=0.82; tox=0.96;
%         fromy=0.22;toy=0.27;
%         uicontrol( ...
%             'Style','push', ...
%             'Units','normalized', ...
%             'BackgroundColor',[0.80 0.80 0.80], ...
%             'ListboxTop',0, ...
%             'Position',[fromx fromy tox-fromx toy-fromy], ...
%             'CallBack',strcat(thisfile,'(''simula'')'), ...
%             'String','Simulate');


        % Compute bounds

%         fromx=0.82; tox=0.96;
%         fromy=0.16;toy=0.21;
%         uicontrol( ...
%             'Style','push', ...
%             'Units','normalized', ...
%             'BackgroundColor',[0.80 0.80 0.80], ...
%             'ListboxTop',0, ...
%             'Position',[fromx fromy tox-fromx toy-fromy], ...
%             'CallBack',strcat(thisfile,'(''bounds'')'), ...
%             'String','Compute Bounds');

        % Compute pt-semiflows

        fromx=0.82; tox=0.96;
        fromy=0.10;toy=0.15;
%         uicontrol( ...
%             'Style','push', ...
%             'Units','normalized', ...
%             'BackgroundColor',[0.80 0.80 0.80], ...
%             'ListboxTop',0, ...
%             'Position',[fromx fromy tox-fromx toy-fromy], ...
%             'CallBack',strcat(thisfile,'(''ptsemis'')'), ...
%             'Tag','ptsemi',...
%             'String','P T semiflows');
        uicontrol( ...
            'Style','push', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'CallBack',strcat(thisfile,'(''simula'')'), ...
            'String','Simulate');

        % Close

        fromx=0.82; tox=0.96;
        fromy=0.04;toy=0.09;
        uicontrol( ...
            'Style','push', ...
            'Units','normalized', ...
            'BackgroundColor',[0.80 0.80 0.80], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'CallBack','close(gcf)', ...
            'String','Close');

        % Bottom frame
        fromx=0.02; tox=0.78;
        fromy=0.02;toy=0.35;
        uicontrol( ...
            'Style','frame', ...
            'Units','normalized', ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'BackgroundColor',[0.70 0.70 0.70]);

        fromx=0.03; tox=0.10;
        fromy=0.3;toy=0.33;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.70 0.70 0.70], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','File:');

        fromx=0.11; tox=0.21;
        fromy=0.29;toy=0.33;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Tag', 'rdpfile', ...
            'String','','Enable','off');

%         fromx=0.22; tox=0.29;
%         fromy=0.3;toy=0.33;
%         uicontrol( ...
%             'Style','push', ...
%             'Units','normalized', ...
%             'BackgroundColor',[0.8 0.8 0.8], ...
%             'ListboxTop',0, ...
%             'Position',[fromx fromy tox-fromx toy-fromy], ...
%             'CallBack',strcat(thisfile,'(''load'')'), ...
%             'Tag','loadFile',...
%             'String','Load');

        fromx=0.23; tox=0.3;
        fromy=0.3;toy=0.33;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.70 0.70 0.70], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','Final Time:');

        fromx=0.3; tox=0.35;
        fromy=0.29;toy=0.33;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Tag', 'tf', ...
            'String','15');

        fromx=0.35; tox=0.47;
        fromy=0.3;toy=0.33;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.70 0.70 0.70], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','Num. simulations:');

        fromx=0.47; tox=0.51;
        fromy=0.29;toy=0.33;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Tag', 'nsim', ...
            'String','200');

        fromx=0.52; tox=0.63;
        fromy=0.30;toy=0.33;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.70 0.70 0.70], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','Hybrid Sampling:');

        fromx=0.63; tox=0.77;
        fromy=0.29;toy=0.33;
        uicontrol( ...
            'Style','popupmenu', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'Value', 2,...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'Tag', 'HybridSampling', ...
            'String','fixed step|variable step'); % |variable step saving data');

        fromx=0.03; tox=0.10;
        fromy=0.25;toy=0.28;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.70 0.70 0.70], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','Pre:');

        fromx=0.11; tox=0.77;
        fromy=0.24;toy=0.28;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'HorizontalAlignment','left', ...
            'Tag', 'Pre', ...
            'String',' [1 0;0 1]',...
            'Callback',strcat(thisfile,'(''import_pre'')'));

        fromx=0.03; tox=0.10;
        fromy=0.2;toy=0.23;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.70 0.70 0.70], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','Post:');

        fromx=0.11; tox=0.77;
        fromy=0.19;toy=0.23;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'HorizontalAlignment','left', ...
            'Tag', 'Post', ...
            'String',' [0 1;1 0]',...
            'Callback',strcat(thisfile,'(''import_post'')'));

        fromx=0.03; tox=0.1;
        fromy=0.15;toy=0.18;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.70 0.70 0.70], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','Mo:');

        fromx=0.11; tox=0.65;
        fromy=0.14;toy=0.18;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'HorizontalAlignment','left', ...
            'Tag', 'M0', ...
            'String',' [2; 1]',...
            'Callback',strcat(thisfile,'(''import_m0'')'));

        fromx=0.03; tox=0.1;
        fromy=0.10;toy=0.13;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.70 0.70 0.70], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','Lambda:');

        fromx=0.11; tox=0.65;
        fromy=0.09;toy=0.13;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'HorizontalAlignment','left', ...
            'Tag', 'Lambda', ...
            'String',' [1; 1]',...
            'Callback',strcat(thisfile,'(''import_lambda'')'));

        fromx=0.03; tox=0.1;
        fromy=0.05;toy=0.08;
        uicontrol( ...
            'Style','text', ...
            'Units','normalized', ...
            'BackgroundColor',[0.70 0.70 0.70], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'String','T. Type:');

        fromx=0.11; tox=0.65;
        fromy=0.04;toy=0.08;
        uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[fromx fromy tox-fromx toy-fromy], ...
            'HorizontalAlignment','left', ...
            'Tag', 'tr_type', ...
            'String',' [''c''; ''c'']',...
            'Callback',strcat(thisfile,'(''import_type'')'));


        set(figpri, 'Visible','on');
        SimHPN('ini_UserData');

    case 'ini_UserData'  % Initialize figure and controls

        data.t = []; % time vector after a simulation
        data.m = []; % marking evolution after a simulation
        data.f = []; % flow evolution after a simulation
        data.m_c = []; % marking evolution applying a control law
        data.w_c = []; % controlled flow evolution corresponding to a control law
        data.u_c = []; % control corresponding to a control law
        data.p_semi = []; % contains the p-semiflows
        data.t_semi = []; % contains the t-semiflows
        data.bounds{1} = [];
        data.bounds{2} = [];
        data.bounds{3} = [];
        pre = eval(get(findobj(gcf,'Tag','Pre'),'String'));
        data.view_m = ones(1,size(pre,1)); % components of m that are represented on the plot
        data.view_f = ones(1,size(pre,2)); % components of f that are represented on the plot
        set(gcf,'UserData',data);

        
        %================
        % LOAD .RDP+.TN + .mat FILE
        %================
    case 'load_rdp'

        [file, path1] = uigetfile({'*.rdp'}, 'Load');
        file = char(file);
        path1 = char(path1);
        file2=fullfile(path1,file);

        if ~isempty(strfind(file,'.rdp'))
            rdpfile=file2;%get(findobj(gcf,'Tag','rdpfile'),'String');
            [Pre,Post,M0, Lambda] = SimHPN_rdp(rdpfile);
            set(findobj(gcf,'Tag','Pre'),'String',sprintf(' %s',mat2str(Pre)));
            set(findobj(gcf,'Tag','Post'),'String',sprintf(' %s',mat2str(Post)));
            set(findobj(gcf,'Tag','M0'),'String',sprintf(' %s',mat2str(M0)));
            set(findobj(gcf,'Tag','Lambda'),'String',sprintf(' %s',mat2str(Lambda)));
            temp = '[';
            for i = 1 : size(Pre,2)-1
                temp = sprintf('%s''c'';',temp);
            end
            temp = sprintf('%s''c'']',temp);
            set(findobj(gcf,'Tag','tr_type'),'String',sprintf(' %s',temp));
            set(findobj(gcf,'Tag','rdpfile'),'String',file);
        end
        SimHPN('ini_UserData');

    case 'load_tn'

        [file, path1] = uigetfile({'*.TN'}, 'Load');
        file = char(file);
        path1 = char(path1);
        file2=fullfile(path1,file);

        if ~isempty(strfind(file,'.TN'))
            [Pre,Post,M0,Lambda]=SimHPN_readTN(file2);
            set(findobj(gcf,'Tag','Pre'),'String',mat2str(Pre));
            set(findobj(gcf,'Tag','Post'),'String',mat2str(Post));
            set(findobj(gcf,'Tag','M0'),'String',mat2str(M0));
            set(findobj(gcf,'Tag','Lambda'),'String',mat2str(Lambda));
            set(findobj(gcf,'Tag','rdpfile'),'String',file);
            temp = '[';
            for i = 1 : size(Pre,2)-1
                temp = sprintf('%s''c'';',temp);
            end
            temp = sprintf('%s''c'']',temp);
            set(findobj(gcf,'Tag','tr_type'),'String',sprintf(' %s',temp));
        end
        SimHPN('ini_UserData');
    case 'load_mat'

        [file, path1] = uigetfile({'*.mat'}, 'Load');
        file = char(file);
        path1 = char(path1);
        file2=fullfile(path1,file);

        if ~isempty(strfind(file,'.mat'))
            eval(sprintf('load %s Pre Post M0 Lambda Trans_Type',file2));
            if (~exist('Pre','var') || ~exist('Post','var') || ~exist('M0','var'))
                h = errordlg('Pre, Post and M0 matrices should be defined in the .mat file','SimHPN Toolbox');
                uiwait(h);
                return;
            end
            set(findobj(gcf,'Tag','Pre'),'String',mat2str(Pre));
            set(findobj(gcf,'Tag','Post'),'String',mat2str(Post));
            set(findobj(gcf,'Tag','M0'),'String',mat2str(M0));
            if exist('Lambda','var')
                set(findobj(gcf,'Tag','Lambda'),'String',mat2str(Lambda));
            else
                temp = '[';
                for i = 1 : size(Pre,2)
                    temp = sprintf('%s''1'';',temp );
                end
                temp = sprintf('%s ]',temp);
                set(findobj(gcf,'Tag','Lambda'),'String',sprintf(' %s',temp));
            end
            if exist('Trans_Type','var')
                temp = '[';
                for i = 1 : size(Pre,2)
                    temp = sprintf('%s''%s'';',temp, Trans_Type(i));
                end
                temp = sprintf('%s ]',temp);
                set(findobj(gcf,'Tag','tr_type'),'String',sprintf(' %s',temp));
            else
                temp = '[';
                for i = 1 : size(Pre,2)
                    temp = sprintf('%s''c'';',temp);
                end
                temp = sprintf('%s ]',temp);
                set(findobj(gcf,'Tag','tr_type'),'String',sprintf(' %s',temp));
            end
            set(findobj(gcf,'Tag','rdpfile'),'String',file);
        end
        SimHPN('ini_UserData');

    case 'save_mat'

        [filename, pathname] = uiputfile('*.mat', 'Save Workspace as');

        Pre=eval(get(findobj(gcf,'Tag','Pre'),'String'));
        Post=eval(get(findobj(gcf,'Tag','Post'),'String'));
        M0=eval(get(findobj(gcf,'Tag','M0'),'String'));
        Lambda=eval(get(findobj(gcf,'Tag','Lambda'),'String'));
        Trans_Type=eval(get(findobj(gcf,'Tag','tr_type'),'String'));

        if (~isequal(filename,0) && ~isequal(pathname,0))
            eval(sprintf('save %s Pre Post M0 Lambda Trans_Type',fullfile(pathname, filename)));
        end
        %================
        % FIRING SEMANTICS
        %================

    case 'infser'

        set(findobj(gcf,'Tag','infser'),'Value',1);
        set(findobj(gcf,'Tag','prod'),'Value',0);

    case 'prod'

        set(findobj(gcf,'Tag','infser'),'Value',0);
        set(findobj(gcf,'Tag','prod'),'Value',1);

        %================
        % PLOT
        %================

    case 'plot_marking'

        set(findobj(gcf,'Tag','pmarking'),'Value',1);
        set(findobj(gcf,'Tag','pflow'),'Value',0);
        set(findobj(gcf,'Tag','pplaces'),'Value',0);
        SimHPN('plot');

    case 'plot_flow'

        set(findobj(gcf,'Tag','pmarking'),'Value',0);
        set(findobj(gcf,'Tag','pflow'),'Value',1);
        set(findobj(gcf,'Tag','pplaces'),'Value',0);
        SimHPN('plot');

    case 'plot_mxm'

        set(findobj(gcf,'Tag','pmarking'),'Value',0);
        set(findobj(gcf,'Tag','pflow'),'Value',0);
        set(findobj(gcf,'Tag','pplaces'),'Value',1);
        SimHPN('plot');

    case 'plot'
        
        data = get(gcf,'UserData');
        t = data.t;
        if isempty(t)
            return;
        end
        m = data.m;
        f = data.f;
        nps = size(m,2);
        nts = size(f,2);
        leg = []; %used for legend command
        %%%%%% check if the selected variables to be plot are consistent
        %%%%%% with the actual model
        if (length(data.view_m) > nps)
            data.view_m = data.view_m(1:nps);
            set(gcf,'UserData',data);
        elseif (length(data.view_m) < nps)
            data.view_m = [data.view_m zeros(1,nps-length(data.view_m))];
            set(gcf,'UserData',data);
        end
        if (length(data.view_f) > nts)
            data.view_f = data.view_f(1:nts);
            set(gcf,'UserData',data);
        elseif (length(data.view_f) < nts)
            data.view_f = [data.view_f zeros(1,nts-length(data.view_f))];
            set(gcf,'UserData',data);
        end
        if get(findobj(gcf,'Tag','pmarking'),'Value')==1
            if isempty(m)
                return;
            end
            drawing = find(data.view_m>0);
            if isempty(drawing)
                warndlg('No place marking is selected to be plot!','SimHPN Toolbox - Warning message');
                return;
            end
            plot(t,m(:,drawing));
            title('Marking Evolution');
            leg(1:length(drawing),1)='m';
            legend([leg int2str((drawing)')]);
        end

        if get(findobj(gcf,'Tag','pflow'),'Value')==1
            if isempty(f)
                return;
            end
            drawing = find(data.view_f>0);
            if isempty(drawing)
                warndlg('No transition flow is selected to be plot!','SimHPN Toolbox - Warning message');
                return;
            end
            plot(t,f(:,drawing));
            title('Flow Evolution');
            leg(1:length(drawing),1)='f';
            legend([leg int2str((drawing)')]);
        end

        if get(findobj(gcf,'Tag','pplaces'),'Value')==1
            if isempty(m)
                return;
            end
            px=eval(get(findobj(gcf,'Tag','px'),'String'));
            py=eval(get(findobj(gcf,'Tag','py'),'String'));
            plot(m(:,px),m(:,py));
            title('Marking vs. Marking');
            xlabel(strcat('m(',num2str(px),')'));
            ylabel(strcat('m(',num2str(py),')'));
        end

        %================
        % Toolbar Visualization
        %================

    case 'opttool'

        if strcmp(get(gcbo, 'Checked'),'on')
            set(gcbo, 'Checked', 'off');
            set(gcf,'ToolBar','none')
        else
            set(gcbo, 'Checked', 'on');
            set(gcf,'ToolBar','figure');
        end

        %================
        % SIMULATE
        %================

    case 'simula'

        Pre=eval(get(findobj(gcf,'Tag','Pre'),'String'));
        Post=eval(get(findobj(gcf,'Tag','Post'),'String'));
        m0=eval(get(findobj(gcf,'Tag','M0'),'String'));
        tf=eval(get(findobj(gcf,'Tag','tf'),'String'));
        lambda=eval(get(findobj(gcf,'Tag','Lambda'),'String'));
        tr_type=eval(get(findobj(gcf,'Tag','tr_type'),'String'));
        nsim=eval(get(findobj(gcf,'Tag','nsim'),'String'));
        mode=get(findobj(gcf,'Tag','HybridSampling'),'Value');
        sampling=eval(get(findobj(gcf,'Tag','sampling'),'String'));
        
        error = 1;
        %Checking for errors
        if (size(Post,1)~=size(Pre,1))||(size(Post,2)~=size(Pre,2))
            errordlg('Matrices Post and Pre do not have the same dimensions','SimHPN: Data error');
        elseif (size(lambda,2)~=1)||(size(lambda,1)~=size(Pre,2))
            errordlg('Vector lambda does not have the correct dimension','SimHPN: Data error');
        elseif (size(tr_type,2)~=1)||(size(tr_type,1)~=size(Pre,2))
            errordlg('Vector type does not have the correct dimension','SimHPN: Data error');
        elseif (size(m0,2)~=1)||(size(m0,1)~=size(Pre,1))
            errordlg('Vector m0 does not have the correct dimension','SimHPN: Data error');
        elseif min(min(Pre))<0
            errordlg('Matrix Pre has negative values','SimHPN: Data error');
        elseif min(m0)<0
            errordlg('Vector m0 has negative values','SimHPN: Data error');
        elseif length(find((tr_type=='q')|(tr_type=='d')|(tr_type=='c')))<length(tr_type)
            errordlg('Type of some transitions wrongly defined','SimHPN: Data error')
        elseif min(lambda)<=0
            errordlg('Vector lambda has non-positive values','SimHPN: Data error');
        elseif min(max(Pre',[],2))==0
            errordlg('According to Matrix Pre, some transition does not have input places','SimHPN: Data error');
        else
            error = 0;
        end;

        if (error == 1)
            return;
        end
        if get(findobj(gcf,'Tag','infser'),'Value')==1
            seman=1;
        else  % Product semantics
            seman=3;
        end

        eabs=eval(get(findobj(gcf,'Tag','eabs'),'String'));
        erel=eval(get(findobj(gcf,'Tag','erel'),'String'));

        %if all transitions are continuous simulate using the old procedure

        tic;
        if (length(find(tr_type=='c'))==length(lambda))
            [m,f,t] = SimHPN_Simul(-21,Pre, Post, m0, tf, lambda, seman, erel, eabs);
        else
            [m,f,t]=SimHPN_SimHyb(-21,Pre,Post,lambda,tr_type,m0,nsim,tf,seman,mode,sampling);
        end
        toc;

        disp(sprintf('Final marking: %s',mat2str(m(size(m,1),:),6)));
        disp(sprintf('Final flow: %s',mat2str(f(size(f,1),:),6)));

        data = get(gcf,'UserData');
        data.t = t;
        data.m = m;
        data.f = f;
        set(gcf,'UserData',data);
        SimHPN('plot');

        %================
        % Compute Bounds
        %================

    case 'bounds'

        Pre=eval(get(findobj(gcf,'Tag','Pre'),'String'));
        Post=eval(get(findobj(gcf,'Tag','Post'),'String'));
        m0=eval(get(findobj(gcf,'Tag','M0'),'String'));
        tf=eval(get(findobj(gcf,'Tag','tf'),'String'));
        lambda=eval(get(findobj(gcf,'Tag','Lambda'),'String'));

        [is_bound_min,is_bound_max, fs_bound] = SimHPN_Bounds ( Pre, Post, m0, lambda, 0);

        data = get(gcf,'UserData');
        data.bounds{1} = is_bound_min;
        data.bounds{2} = is_bound_max;
        data.bounds{3} = fs_bound;
        set(gcf,'UserData',data);
        
        str = sprintf('Throughput lower bound under Infinite server semantics:\n\n %s\n\n Throughput upper bound under Infinite server semantics: \n \n %s \n\n Throughput bound under Finite server semantics: %s \n\n', ...
            mat2str(is_bound_min), mat2str(is_bound_max), mat2str(fs_bound));
        h = figure('WindowStyle','modal','Units','normalized','Menubar','None','Name','SimHPN Toolbox',...
            'Position',[0.35 0.4 0.5 0.4],...
            'NumberTitle','off','CloseRequestFcn','delete(gcf)');
        uicontrol('Units','normalized','Style','text','Position',[0 0.15 1 0.85],'HorizontalAlignment','left',...
            'String',str);
        uicontrol('Units','normalized','Style','pushbutton','String','Close','Callback','delete(gcf)',...
            'Position',[0.2 0.02 0.2 0.1],'HorizontalAlignment','left');
        uicontrol('Units','normalized','Style','pushbutton','String','Save to Workspace','Callback',strcat(thisfile,'(''save_bounds'')'),...
            'Position',[0.6 0.02 0.2 0.1],'HorizontalAlignment','left');
    case 'save_bounds'
        delete(gcf);
        data = get(gcf,'UserData');
        checkLabels = {'Save the throughput lower bound for Infinite server semantics:', ...
            'Save the throughput upper bound for Infinite server semantics:',...
            'Save the throughput bound for finite server semantics:'};
        varNames = {'is_bound_min', 'is_bound_max','fs_bound'};
        items = {data.bounds{1}, data.bounds{2}, data.bounds{3}};
        export2wsdlg(checkLabels, varNames, items, 'Save the bounds to Workspace');

        
        %================
        % Check conservativeness
        %================
    case 'conservativeness'
        Pre = eval(get(findobj(gcf,'Tag','Pre'),'String'));
        Post = eval(get(findobj(gcf,'Tag','Post'),'String'));
        % create the constraints matrix
        a = (Post-Pre)';
        b = zeros(size(Pre,2),1);
        ctype = '';
        for i = 1 : size(Pre,2)
            ctype = sprintf('%sS',ctype);;
        end
        vartype = '';
        vartype2 = '';
        for i = 1 : size(Pre,1)
            vartype = sprintf('%sI',vartype);;
            vartype2 = sprintf('%sC',vartype2);;
        end
        lb = ones(size(Pre,1),1);
        c = lb;
        [xmin, fmin, status, extra] = glpk (c, a, b, lb, [], ctype, vartype);
        
        if (status == 5)
            uiwait(msgbox(sprintf('The net is conservative!\n\n A conservative component: %s\n (see the MATLAB command prompt)',mat2str(xmin)),'SimHPN Toolbox','modal'));
            fprintf(1,'The net is conservative. A conservative component:\n%s\n',mat2str(xmin));
        else
            ButtonName = questdlg('The net is not conservative.\nCompute the places not belonging to any P-semiflow support?', ...
                'SimHPN Toolbox', 'Yes', 'No', 'No');
            if strcmpi(ButtonName,'No')
                return;
            else
                str = 'Places not belonging to any P-semiflow support:';
                for i = 1 : size(Pre,1)
                    lb = zeros(size(Pre,1),1);
                    lb(i) = 1;
                    [xmin, fmin, status, extra] = glpk (c, a, b, lb, [], ctype, vartype2);
                    if (status ~= 5)
                        str = sprintf('%s p_%d',str,i);
                    end
                end
                uiwait(msgbox(str,'SimHPN Toolbox','modal'));
                fprintf(1,'\n%s\n',str);
            end
        end        
        
        %================
        % Check consistency
        %================
    case 'consistency'
        Pre = eval(get(findobj(gcf,'Tag','Pre'),'String'));
        Post = eval(get(findobj(gcf,'Tag','Post'),'String'));
        % create the constraints matrix
        a = (Post-Pre);
        b = zeros(size(Pre,1),1);
        ctype = '';
        for i = 1 : size(Pre,1)
            ctype = sprintf('%sS',ctype);;
        end
        vartype = '';
        vartype2 = '';
        for i = 1 : size(Pre,2)
            vartype = sprintf('%sI',vartype);;
            vartype2 = sprintf('%sC',vartype2);;
        end
        lb = ones(size(Pre,2),1);
        c = lb;
        [xmin, fmin, status, extra] = glpk (c, a, b, lb, [], ctype, vartype);
        
        if (status == 5)
            uiwait(msgbox(sprintf('The net is consistent!\n\n A repetitive component: %s\n (see the MATLAB command prompt)',mat2str(xmin)),'SimHPN Toolbox','modal'));
            fprintf(1,'The net is consistent. A repetitive component:\n%s\n',mat2str(xmin));
        else
            ButtonName = questdlg('The net is not consistent.\nCompute the transitions not belonging to any T-semiflow support?', ...
                'SimHPN Toolbox', 'Yes', 'No', 'No');
            if strcmpi(ButtonName,'No')
                return;
            else
                str = 'Transitions not belonging to any T-semiflow support:';
                for i = 1 : size(Pre,2)
                    lb = zeros(size(Pre,2),1);
                    lb(i) = 1;
                    [xmin, fmin, status, extra] = glpk (c, a, b, lb, [], ctype, vartype2);
                    if (status ~= 5)
                        str = sprintf('%s t_%d',str,i);
                    end
                end
                uiwait(msgbox(str,'SimHPN Toolbox','modal'));
                fprintf(1,'\n%s\n',str);
            end
        end        

        
        %================
        % P-T Semiflows
        %================

    case 'ptsemis'

        Pre=eval(get(findobj(gcf,'Tag','Pre'),'String'));
        Post=eval(get(findobj(gcf,'Tag','Post'),'String'));
        [p_semi,t_semi] = SimHPN_ptsemis(Post-Pre);
        data = get(gcf,'UserData');
        data.p_semi = p_semi;
        data.t_semi = t_semi;
        set(gcf,'UserData',data);

        str = sprintf('P-Semiflows (each column is a P-semiflow):\n\n %s \n\n T-Semiflows (each column is a T-semiflow):\n\n %s\n',mat2str(p_semi),mat2str(t_semi));
        h = figure('WindowStyle','modal','Units','normalized','Menubar','None','Name','SimHPN Toolbox',...
            'Position',[0.35 0.4 0.5 0.4],...
            'NumberTitle','off','CloseRequestFcn','delete(gcf)');
        uicontrol('Units','normalized','Style','text','Position',[0 0.15 1 0.85],'HorizontalAlignment','left',...
            'String',str);
        uicontrol('Units','normalized','Style','pushbutton','String','Close','Callback','delete(gcf)',...
            'Position',[0.2 0.02 0.2 0.1],'HorizontalAlignment','left');
        uicontrol('Units','normalized','Style','pushbutton','String','Save to Workspace','Callback',strcat(thisfile,'(''save_ptsemi'')'),...
            'Position',[0.6 0.02 0.2 0.1],'HorizontalAlignment','left');
    case 'save_ptsemi'
        delete(gcf);
        data = get(gcf,'UserData');
        checkLabels = {'Save the P-semiflows to variable named:' ...
            'Save the T-semiflows to variable named:'};
        varNames = {'p_semi', 't_semi'};
        items = {data.p_semi, data.t_semi};
        export2wsdlg(checkLabels, varNames, items, 'Save P T Semiflows to Workspace');

        %========================
        % Compute control laws %
        %========================
    case 'contr_law'
        Pre = eval(get(findobj(gcf,'Tag','Pre'),'String'));
        Post = eval(get(findobj(gcf,'Tag','Post'),'String'));
        lambda = eval(get(findobj(gcf,'Tag','Lambda'),'String'));
        m0 = eval(get(findobj(gcf,'Tag','M0'),'String'));
        sampling = eval(get(findobj(gcf,'Tag','sampling'),'String'));
        
        [ s, m, u, w ] = SimHPN_control_law( Pre , Post , lambda , m0 , sampling, inf );

         data = get(gcf,'UserData');
         data.m_c = m;
         data.u_c = u;
         data.w_c = w;
         set(gcf,'UserData',data);
 
         disp(sprintf('Final marking is reached in %d steps, using %s method',s, inf));        
    case 'save_control'
        data = get(gcf,'UserData');
        if isempty(data.m_c)
            errordlg('Compute a control law first!');
            return;
        end
        checkLabels = {'Save the marking evolution to variable named:' ...
            'Save the controlled flow evolution to variable named:'...
            'Save the control input to variable named:'};
        varNames = {'m_c', 'w_c','u_c'};
        items = {data.m_c, data.w_c, data.u_c};
        export2wsdlg(checkLabels, varNames, items, 'Save Control Law Computation to Workspace');
        
        %========================
        % Optimal Observability %
        %========================
    case 'optobs'
        Pre=eval(get(findobj(gcf,'Tag','Pre'),'String'));
        Post=eval(get(findobj(gcf,'Tag','Post'),'String'));
        lambda=eval(get(findobj(gcf,'Tag','Lambda'),'String'));
        temp = mat2str(ones(1,size(Pre,1)));
        cost =char(inputdlg('Measuring cost of places:','Measuring cost of places',1,{temp}));
        if isempty(cost)
            return;
        end
        cost = eval(cost);
        SimHPN_OptObs(Pre,Post,lambda,cost);
        %========================
        % Influence of the controllable transitions %
        %========================
    case 'influ'
        Pre=eval(get(findobj(gcf,'Tag','Pre'),'String'));
        Post=eval(get(findobj(gcf,'Tag','Post'),'String'));
        temp = mat2str([]);
        Tc =char(inputdlg('Set of controllable transitions:','Tc',1,{temp}));
%         if isempty(Tc)
%             return;
%         end
        Tc = eval(Tc);
        SimHPN_Influ(Pre,Post,Tc);
        %========================
        % Net rank-controllability test %
        %========================
    case 'nrc'
        Pre=eval(get(findobj(gcf,'Tag','Pre'),'String'));
        Post=eval(get(findobj(gcf,'Tag','Post'),'String'));
        temp = mat2str([]);
        Tc =char(inputdlg('Set of controllable transitions:','Tc',1,{temp}));
%         if isempty(Tc)
%             return;
%         end
        Tc = eval(Tc);
        SimHPN_NRC(Pre,Post,Tc);
        %================
        % Optimal steady-state
        %================
    case 'optsteady'

        Pre=eval(get(findobj(gcf,'Tag','Pre'),'String'));
        Post=eval(get(findobj(gcf,'Tag','Post'),'String'));
        m0=eval(get(findobj(gcf,'Tag','M0'),'String'));
        tf=eval(get(findobj(gcf,'Tag','tf'),'String'));
        lambda=eval(get(findobj(gcf,'Tag','Lambda'),'String'));

        [is_bound_min,is_bound_max, fs_bound] = SimHPN_Bounds ( Pre, Post, m0, lambda, 1);
        if (isempty(is_bound_min) && isempty(is_bound_max) && isempty(fs_bound))
            return;
        end
        nps = size(Pre,1);
        nts = size(Pre,2);
        disp('Optimal Steady-State under Infinite server semantics:');
        disp(is_bound_max);
        disp('Optimal marking:');
        disp(is_bound_min(nts+1:nts+nps));
        disp('Control input:');
        disp(is_bound_min(nts+nps+1:2*nts+nps));
        
        %================
        % Diagnosis of untimed continuous Pteri net
        %================
    case 'diagnosis'
        Pre=eval(get(findobj(gcf,'Tag','Pre'),'String'));
        Post=eval(get(findobj(gcf,'Tag','Post'),'String'));
        M0=eval(get(findobj(gcf,'Tag','M0'),'String'));

        answer = inputdlg({...
            sprintf('Number of observable transitions (first xxx transitions will be observable):'),...
            sprintf('Sequence of observed transitions:'),...
            sprintf('Sequence of observed firing quantity of transitions:')...
            sprintf('Number of fault classes: \n\t')},'SimHPN: Diagnosis using Continous Petri nets',...
            [1,1,1,1],{'1','[1]','[1]','1'});
        try
            obsT = eval(answer{1});
            obsW = eval(answer{2});
            obsA = eval(answer{3});
            nFaults = eval(answer{4});
        catch except
            errordlg('Error introducing datas!','SimHPN Toolbox');
            return;
        end
        if isempty(nFaults)
            errordlg('Please introduce the number of falt classes.');
            return;
        end
        Tf = {};
        i = 1;
        while (i <= nFaults)
            answer = inputdlg(sprintf('Transitions in fault class %d: \n \t',i),'SimHPN toolbox',1,{'[]'});
            if isempty(answer)
                return;
            end
            if isempty(eval(answer{1}))
                h = errordlg('Please introduce the transitions belonging to the fault class.','SimHPN Toolbox');
                uiwait(h);
            else
                Tf{i} = eval(answer{1});
                if ((min(Tf{i}) <= obsT) || (max(Tf{i}) > size(Pre,2)))
                    h = errordlg(sprintf('t%d is not an unobservable transition!',Tf{i}),'SimHPN Toolbox');
                    uiwait(h);
                else
                    i = i + 1;
                end
            end
        end
        fprintf(1,'=========================================\n');
        fprintf(1,'\nNumber of observable transitions: %d \nFirst %d transitions are observable and the rest not\n',obsT,obsT);
        
        for i = 1 : length(Tf)
            fault = Tf{i};
            str = sprintf('Fault class Tf^{%d}={\\epsilon_',i);
            for j = 1 : length(fault) - 1
                str = sprintf('%s%d, \\epsilon_',str,fault(j));
            end
            str = sprintf('%s%d}',str,fault(length(fault)));
            fprintf(1,'%s\n',str);
        end
        
        if isempty(obsW) %compute only the consistent marking for empty word and return
            SimHPN_diagnoser(Pre,Post,M0,obsT,obsW,obsA,Tf);
            return
        end
        
        if (length(obsW) ~= length(obsA))
            error('The size of obseved transition vector and the observed quantity should be the same');
        end
        
        str='t';
        for i = 1 : length(obsW)-1
            str = sprintf('%s%d(%s)t',str,obsW(i),num2str(obsA(i)));
        end
        str = sprintf('%s%d(%s)',str,obsW(length(obsW)),num2str(obsA(length(obsA))));
        
        fprintf(1,'\nObserved sequence: %s \n',str);
        
        fprintf(1,'=========================================\n');
        fprintf(1,'==========press a key to continue========\n');
        fprintf(1,'=========================================\n\n');
        
        pause
        SimHPN_diagnoser(Pre,Post,M0,obsT,obsW,obsA,Tf);

        %========================
        % MENU SIMULATION %
        %========================
        
    case 'save_simul'
        
        data = get(gcf,'UserData');
        if isempty(data.t)
            errordlg('Do a simulation first!');
            return;
        end
        checkLabels = {'Save time vector to variable named:' ...
            'Save marking evolution matrix to variable named:'...
            'Save flow evolution matrix to variable named:'};
        varNames = {'t', 'm','f'};
        items = {data.t, data.m, data.f};
        export2wsdlg(checkLabels, varNames, items, 'Save Simulation Results to Workspace');
    case 'plot_var_m'
        data = get(gcf,'UserData');
        d = find(data.view_m > 0);
        Pre=eval(get(findobj(gcf,'Tag','Pre'),'String'));
        nps = size(Pre,1);
        str = {};
        for i = 1 : nps
            str{i} = sprintf('m %d',i);
        end
        [s,ok] = listdlg('PromptString','Select the markings to be plot:',...
            'SelectionMode','multiple','InitialValue',d,'Name','Select the desired places',...
            'ListString',str);
        if (ok == 0)
            return;
        end
        new_draw = zeros(1,nps);
        new_draw(s) = 1;
        data.view_m = new_draw;
        set(gcf,'UserData',data);
    case 'plot_var_f'
        data = get(gcf,'UserData');
        d = find(data.view_f > 0);
        Pre=eval(get(findobj(gcf,'Tag','Pre'),'String'));
        nts = size(Pre,2);
        str = {};
        for i = 1 : nts
            str{i} = sprintf('f %d',i);
        end
        [s,ok] = listdlg('PromptString','Select the flows to be plot:',...
            'SelectionMode','multiple','InitialValue',d,'Name','Select the desired transitions',...
            'ListString',str);
        if (ok == 0)
            return;
        end
        new_draw = zeros(1,nts);
        new_draw(s) = 1;
        data.view_f = new_draw;
        set(gcf,'UserData',data);
        
        
    case 'import_pre' %import pre matrix from a workspace defined variable
        str = get(findobj(gcf,'Tag','Pre'),'String');
        Pre=evalin('caller',str);
        if isnumeric(Pre)
            set(findobj(gcf,'Tag','Pre'),'String',mat2str(Pre));
            return;
        end
        Pre=evalin('base',str);
        if isnumeric(Pre)
            set(findobj(gcf,'Tag','Pre'),'String',mat2str(Pre));
        end
        SimHPN('ini_UserData');

    case 'import_post' %import post matrix from a workspace defined variable
        str = get(findobj(gcf,'Tag','Post'),'String');
        Post=evalin('caller',str);
        if isnumeric(Post)
            set(findobj(gcf,'Tag','Post'),'String',mat2str(Post));
            return;
        end
        Post=evalin('base',str);
        if isnumeric(Post)
            set(findobj(gcf,'Tag','Post'),'String',mat2str(Post));
        end
        SimHPN('ini_UserData');

    case 'import_m0' %import m0 matrix from a workspace defined variable
        str = get(findobj(gcf,'Tag','M0'),'String');
        m0 = evalin('caller',str);
        if isnumeric(m0)
            set(findobj(gcf,'Tag','M0'),'String',mat2str(m0));
            return;
        end
        m0 = evalin('base',str);
        if isnumeric(m0)
            set(findobj(gcf,'Tag','M0'),'String',mat2str(m0));
        end
        SimHPN('ini_UserData');
        
    case 'import_lambda' %import lambda matrix from a workspace defined variable
        str = get(findobj(gcf,'Tag','Lambda'),'String');
        lam = evalin('caller',str);
        if isnumeric(lam)
            set(findobj(gcf,'Tag','Lambda'),'String',mat2str(lam));
            return;
        end
        lam = evalin('base',str);
        if isnumeric(lam)
            set(findobj(gcf,'Tag','Lambda'),'String',mat2str(lam));
        end
        SimHPN('ini_UserData');

    case 'import_type' %import transition types from a workspace defined variable
        str = get(findobj(gcf,'Tag','tr_type'),'String');
        type = evalin('caller',str);
        if isnumeric(type)
            set(findobj(gcf,'Tag','tr_type'),'String',mat2str(type));
            return;
        end
        type = evalin('base',str);
        if isnumeric(type)
            set(findobj(gcf,'Tag','tr_type'),'String',mat2str(type));
        end
        SimHPN('ini_UserData');
    case 'PetriNetToolbox'
        Pre = varargin{1};
        Post = varargin{2};
        M0 = varargin{3};
        Lambda = varargin{4};
        file = varargin{5};
        set(findobj(gcf,'Tag','Pre'),'String',sprintf(' %s',mat2str(Pre)));
        set(findobj(gcf,'Tag','Post'),'String',sprintf(' %s',mat2str(Post)));
        set(findobj(gcf,'Tag','M0'),'String',sprintf(' %s',mat2str(M0)));
        set(findobj(gcf,'Tag','Lambda'),'String',sprintf(' %s',mat2str(Lambda)));
        set(findobj(gcf,'Tag','rdpfile'),'String',file);
        temp = '[';
        for i = 1 : size(Pre,2)-1
            temp = sprintf('%s''c'';',temp);
        end
        temp = sprintf('%s''c'']',temp);
        set(findobj(gcf,'Tag','tr_type'),'String',sprintf(' %s',temp));
        
        
        
end;    % switch
