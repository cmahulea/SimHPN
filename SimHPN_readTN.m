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

function [Pre,Post,M0,Lambda,placesName,transName]=readTN(file2);

try
    a=textread(file2,'%s');     %read .TN files
catch
    disp(lasterr);
    return;
end

nts = [];
nps = [];
for i = 1 : length(a)
    if (~isempty(strfind(lower(a{i}),'places:')) & isempty(nps))
        nps = str2num(a{i+1});
    end
    if (~isempty(strfind(lower(a{i}),'transitions:')) & (isempty(nts)))
        nts = str2num(a{i+1});
        break;
    end
end
a(1:i) = [];
markpar.name={};
markpar.val=[];
delaypar.name={};
delaypar.val=[];
startPlaces = 0;
startTransitions = 0;
places.name={};
places.m0=[];
transitions = [];
lama = [];
for i = 1 : length(a)
    if ~isempty(strfind(lower(a{i}),'markpar'))
        markpar.name{length(markpar.name)+1} = a{i+1};
        markpar.val = [markpar.val str2num(a{i+2})];
    end
    if ~isempty(strfind(lower(a{i}),'delaypar'))
        delaypar.name{length(delaypar.name)+1} = a{i+1};
        delaypar.val = [delaypar.val str2num(a{i+2})];
    end
    
    %read places
    if (~isempty(strfind(lower(a{i}),'place')) & (startPlaces==1))
        places.name{length(places.name)+1} = a{i+1};
        temp = str2num(a{i+2});
        if ~isempty(temp)
            places.m0 = [places.m0 temp];
        else
            for j = 1 : length(markpar.name)
                if strcmp(markpar.name{j},a{i+2})
                    places.m0 = [places.m0 markpar.val(j)];
                    break;
                end
            end
        end
    end
    
    %read transitions
    if (~isempty(strfind(lower(a{i}),'transition')) && (startTransitions==1))
        transitions(length(transitions)+1).name = a{i+1};
        transitions(length(transitions)).lam = [];
        temp = str2num(a{i+2});
        if strcmp(a{i+2},'<MD>')
            temp = 1;
        end
        if ~isempty(temp)
            transitions(length(transitions)).lam = 1/temp;
            lama = [lama transitions(length(transitions)).lam];
        else
            for j = 1 : length(delaypar.name)
                if strcmp(delaypar.name{j},a{i+2})
                    transitions(length(transitions)).lam = 1/delaypar.val(j);
                    lama = [lama transitions(length(transitions)).lam];
                    break;
                end
            end
        end
        if (strcmp(a{i+4},'IM'))
            transitions(length(transitions)).im = 1;
            lama(length(lama)) = Inf;
        else
            transitions(length(transitions)).im = 0;
        end
        transitions(length(transitions)).inArcs={};
        transitions(length(transitions)).outArcs={};
        transitions(length(transitions)).inArcsV=[];
        transitions(length(transitions)).outArcsV=[];
        
        if ~isempty(transitions(length(transitions)).lam)  %read the arcs of the transition
            index = i;
            cont = 1;
            while (cont == 1)
                if (strcmp(a{index},'INPARCS'))
                    ninArcs = str2num(a{index+1});
                    inArcs = index + 2;
                end
                if (strcmp(a{index},'OUTPARCS'))
                    noutArcs = str2num(a{index+1});
                    outArcs = index + 2;
                    cont = 0;
                end
                index = index + 1;
            end
            
            % read input arcs
            index = inArcs;
            for j = 1 : ninArcs
                transitions(length(transitions)).inArcsV = [transitions(length(transitions)).inArcsV str2num(a{index})];
                transitions(length(transitions)).inArcs{j} = a{index+1};
                salt = str2num(a{index+2});
                index = index + 2*salt + 3;
            end
            % read output arcs
            index = outArcs;
            for j = 1 : noutArcs
                transitions(length(transitions)).outArcsV = [transitions(length(transitions)).outArcsV str2num(a{index})];
                transitions(length(transitions)).outArcs{j} = a{index+1};
                salt = str2num(a{index+2});
                index = index + 2*salt + 3;
            end
        end
    end
    
    if (~isempty(strfind(lower(a{i}),'list')) && ~isempty(strfind(lower(a{i+1}),'of')) ...
            && ~isempty(strfind(lower(a{i+2}),'places')))
        startPlaces=1;
    end
    if (~isempty(strfind(lower(a{i}),'list')) && ~isempty(strfind(lower(a{i+1}),'of')) ...
            && ~isempty(strfind(lower(a{i+2}),'transitions')))
        startPlaces=0;
        startTransitions = 1;
    end
end

disp(sprintf('Lambda: %s',mat2str(lama,2)));
places.name(2) = [];
places.name(1) = [];
transitions(2) = [];
transitions(1) = [];
Pre = zeros(nps,nts);
Post = zeros(nps,nts);
M0 = zeros(nps,1);
Lambda = zeros(nts,1);

%print transitions & places name
str = 'Transitions = [';
transName = {};
for i = 1 : nts - 1
    transName{i} = transitions(i).name;
    str = sprintf('%s%s,',str,transitions(i).name);
end
transName{nts} = transitions(nts).name;
str = sprintf('%s%s]',str,transitions(nts).name);
disp(str);

placesName={};
str = 'Places = [';
for i = 1 : nps - 1
    placesName{i} = places.name{i};
    str = sprintf('%s%s,',str,places.name{i});
end
placesName{nps} = places.name{nps};
str = sprintf('%s%s]',str,places.name{nps});
disp(str);

for i = 1 : nts
    if (transitions(i).im == 1)
        Lambda(i) = Inf;
    else
        Lambda(i) = transitions(i).lam;
    end
    for j = 1 : length(transitions(i).inArcs)
        for k = 1 : length(places.name)
            if strcmp(places.name{k},transitions(i).inArcs{j})
                Pre(k,i) = transitions(i).inArcsV(j);
                break;
            end
        end
    end
    for j = 1 : length(transitions(i).outArcs)
        for k = 1 : length(places.name)
            if strcmp(places.name{k},transitions(i).outArcs{j})
                Post(k,i) = transitions(i).outArcsV(j);
                break;
            end
        end
    end
end

disp(sprintf('M0: %s',mat2str(places.m0,2)));
M0 = places.m0';
return
