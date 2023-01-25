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

function [ret,costret] = SimHPN_OptObs(Pre,Post,lambda,cost)

ret = [];
costret = 0;

%cut the synchronizations
remove = [];
for i = 1 : size(Pre,2)
    p = find(Pre(:,i));
    if length(p) > 1
        remove = [remove i];
    end
end
% if there exists a place with input and output transition a join measure
% this place
for i = 1 : size(Pre,1)
    tin = find(Post(i,:));
    tout = find(Pre(i,:));
    if isempty(setdiff(union(tin,tout),remove))
        ret = [ret i];
        costret = costret + cost(i);
    end
end
Pre(:,remove) = [];
Post(:,remove) = [];
lambda(remove) = [];

if (isempty(Pre) && isempty(Post))
    return;
end

%remove the choices
remove = [];
for i = 1 : size(Pre,1)
    t = find(Pre(i,:));
    if (length(t) <= 1)
        continue;
    end
    Post = [Post zeros(size(Post,1),1)];
    Pre = [Pre zeros(size(Pre,1),1)];
    Pre(i,size(Pre,2)) = 1;
    for k = 1 : length(t)
        p = find(Post(:,t(k)));
        for j = 1 : length(p)
            Post(p(j),size(Post,2)) = (Post(p(j),t(k))/Pre(i,t(k))) * (lambda(t(k)) / sum(lambda(t)));
        end
    end
    lambda = [lambda sum(lambda(t))];
    remove = [remove t];
end
Pre(:,remove) = [];
Post(:,remove) = [];
lambda(remove) = [];

% compute F (i.e. terminal places)
F = [];
for i = 1 : size(Pre,1)
    tout = find(Pre(i,:));
    if isempty(tout)
        F = [F i];
        continue;
    end
    pout = [];
    for j = 1 : length(tout)
        pout = [pout find(Post(:,tout(j)))];
    end
    if isempty(pout)
        F = [F i];
    end
end

%compute the attributions
attrib = [];
for i = 1 : size(Post,1)
    temp = find(Post(i,:));
    for j = length(temp) : -1 : 1
        if isempty(find(Pre(:,temp(j))))
            temp(j) = [];
        end
    end
    if length(temp) > 1
        attrib = [attrib i];
    end
end

%if is a cycle return the solution
if (isempty(attrib) && isempty(F))
    [y,i] = min(cost);
    ret = i(1);
    costret = cost(i(1));
    [tt,ii] = find(cost == costret);
    disp(sprintf('Minimum cost is: %s',num2str(costret)));
    for j = 1 : length(tt)
        disp(sprintf('Solution %s:',num2str(j)));
        disp(sprintf('Measured place is: p%s',num2str(ii(j))));
    end
    return;
end

%if is a AF the solution is F
if isempty(attrib)
    ret = union(ret,F);
    costret = costret + sum(cost(F));
    disp(sprintf('Minimum cost is: %s',num2str(costret)));
    disp('Measured places are:');
    disp(ret);
    return;
end

% compute the places estimated from F
GF = SimHPN_Estimates(F,Pre,Post);

%compute the incoming branches
G = {};
for i = 1 : length(attrib)
    tin = find(Post(attrib(i),:));
    for j = 1 : length(tin)
        pin = find(Pre(:,tin(j)));
        temp = SimHPN_Estimates(pin,Pre,Post);
        remove = [];
        for k = length(temp) : -1 : 1
            for l = 1 : length(GF)
                if temp(k) == GF(l)
                    remove = [remove k];
                end
            end
        end
        temp(remove) = [];
        G{length(G)+1} = temp;
    end
end

% update the threads looking at the cost
for i = 1 : length(G)
    temp = G{i};
    remove = [];
    for j = 1 : length(temp)-1
        for k = j+1 : length(temp)
            if (cost(temp(j)) <= cost(temp(k)))
                remove = union(remove,k);
            end
        end
    end
    temp(remove) = [];
    G{i} = temp;
end


%make them disjoints
for i = 1 : length(G)-1
    for j = i+1 : length(G)
        temp1 = G{i};
        temp2 = G{j};
        remove = [];
        for k = length(temp2): -1 : 1
            for l = 1 : length(temp1)
                if (temp2(k) == temp1(l))
                    remove = [remove k];
                end
            end
        end
        temp2(remove) = [];
        G{j} = temp2;
    end
end

%compute A
A = SimHPN_computeA(Pre,Post,lambda);

for i = length(G): -1 : 1
    if isempty(G{i})
        G(i) = [];
    else
        temp = G{i};
        temp = [temp 0];
        G{i} = temp;
    end
end

tic;

[p,c] = SimHPN_Generate(G,A,F,cost);

toc;

costret = c + sum(cost(F));
ret = union(p,F);
disp(sprintf('Minimum cost is: %s',num2str(costret)));
disp('Measured places are:');
disp(ret);

return

