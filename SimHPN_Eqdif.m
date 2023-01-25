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

function dm = SimHPN_Eqdif(t, m, options, Pre, Post, C, lambda, seman, eabs, tf)

persistent h_handle;

if ((nargin == 1) && ishandle(h_handle))
    delete(gcf);
    return;
end
if ((isempty(h_handle)) && (t < tf))
    %   h_handle=waitbar(0,'Solving ODE...');
    h_handle = waitbar(0,'Please wait...','Name','SimHPN Toolbox','CreateCancelBtn','SimHPN_Eqdif(''back'')',...
        'CloseRequestFcn','SimHPN_Eqdif(''back'')');
elseif ~isempty(h_handle)
    %waitbar(t/tf,h_handle,'Solving ODE...');
    waitbar(t/tf,h_handle,'Please wait...');
end

if ((t>=tf) && (~isempty(h_handle)))
    close(h_handle);
    h_handle=[];
end

switch seman
    case 1,    % Infinite servers semantics

        for i=1:size(Pre,2)		% Flow for each transition
            sens = -1;
            for j=1:size(Pre,1)	% Enabling degree
                if ( Pre(j,i) > 0 )
                    if ( sens == -1 ) % First input place
                        sens=m(j)/Pre(j,i);
                    else
                        sens = min(sens, m(j)/Pre(j,i));
                    end
                end
            end
            f(i,1) = lambda(i)*sens;
        end

    case 3,    % Product firing semantics

        for i=1:size(Pre,2)
            sens = -1;
            for j=1:size(Pre,1)
                if ( Pre(j,i) > 0 )
                    if ( sens == -1 )
                        sens=m(j)/Pre(j,i);
                    else
                        sens = sens*m(j)/Pre(j,i);
                    end
                end
            end
            f(i,1) = lambda(i)*sens;
        end

end

dm = C * f;

% Last components of dm are used to store the evolution of
% the throughput of the transitions

dm=[dm;f];
