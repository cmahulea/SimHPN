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

function [m,f,t] = SimHPN_Simul(demo,Pre, Post, m0, tf, lambda, seman, erel, eabs)

if (nargin < 1)
    errordlg('This is a demo version. Please contact us by email at simhpn@unizar.es for the full version.','SimHPN Toolbox for MATLAB');
    return;
end
if (demo ~= -21)
    errordlg('This is a demo version. Please contact us by email at simhpn@unizar.es for the full version.','SimHPN Toolbox for MATLAB');
    return;
end
if (nargin ~= 9)
    errordlg('The number of input parameters is not correct!.','SimHPN Toolbox for MATLAB');
    return;
end

%Checking for errors
if (size(Post,1)~=size(Pre,1))||(size(Post,2)~=size(Pre,2))
    errordlg('Matrices Post and Pre do not have the same dimensions','Data error');  
    return;
elseif (size(lambda,2)~=1)||(size(lambda,1)~=size(Pre,2))
    errordlg('Vector lambda does not have the correct dimension','Data error');
    return;
elseif (size(m0,2)~=1)||(size(m0,1)~=size(Pre,1))
    errordlg('Vector m0 does not have the correct dimension','Data error');
    return;
elseif min(min(Pre))<0
    errordlg('Matrix Pre has negative values','Data error');
    return;
elseif min(min(Post))<0
    errordlg('Matrix Post has negative values','Data error');
    return;
elseif min(m0)<0
    errordlg('Vector m0 has negative values','Data error');
    return;
elseif min(lambda)<=0
    errordlg('Vector lambda has non-positive values','Data error');
    return;
elseif min(max(Pre',[],2))==0
    errordlg('According to matrix Pre, some transition does not have input places','Data error');
    return;
end;

if (nargin < 2)
    error('Use: SimHPN_Simul(Pre, Post,[m0],[tf],[lambda],[seman],[erel],[eabs])');
end

if (seman == 2)
    seman = 3;
end

C = Post - Pre;  % Incidence matrix

if (exist('m0','var') ~= 1) 
    m0 = ones(size(C,1),1); 
end
% if (size(m0,1) == 1) m0 = m0'; end			% Column vector
% if (exist('tf','var') ~= 1) tf = 20; end
% if (exist('lambda','var') ~= 1) lambda = ones(size(C,2),1); end
% if (size(lambda,1) == 1) lambda = lambda'; end % Column vector
% if (exist('seman','var') ~= 1) seman = 1; end
% if (exist('erel','var') ~= 1) erel = 1e-3; end
% if (exist('eabs','var') ~= 1) eabs = 1e-3; end


% Solve differential equation system

nts=size(C,2); % Number of transitions
nps=size(C,1); % Number of places

m0=[m0;zeros(nts,1)];
options = odeset('RelTol',erel,'AbsTol',eabs);
[t,m] = ode45('SimHPN_Eqdif',[0,tf], m0,options, Pre, Post, C, lambda, seman, eabs, tf);



% The first columns of m correspond to the marking of places. The
%  rest are the integral of the throughput of transitions.

% Get flow of transitions

fd=diff(m(:,nps+1:nps+nts));
td=diff(t);

for i=1:nts
    f(:,i)=fd(:,i)./td;
end

% Decouple marking from flow

m=m(:,1:nps);
if ((size(f,1)+1)==size(m,1))
    f = [f(1,:) ; f];
end
