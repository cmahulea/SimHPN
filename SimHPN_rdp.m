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

function [Pre,Post,M0, lambda] = SimHPN_rdp(filename)
% [Pre,Post,M0] = rdp(filename)
% reads matrices Pre, Post and markings from "PM Editeur" format file
% Calls - getnum.m
fid = fopen(filename,'r');
F = fread(fid);
s = setstr(F');
%********** find parameters **********
pla=findstr(s,'NB_PLACE');
pl=getnum(s,pla);
tra=findstr(s,'NB_TRANS');
tr=getnum(s,tra);
arc=findstr(s,'NB_ARC');
ar=getnum(s,arc);
%********** find marking **********
for k=1:pl
    elpl=sprintf('[P%d]',k);
    xpl=findstr(s,elpl);
    mar=findstr(s(xpl:length(s)),'MARKING');
    ma(k)=getnum(s,xpl+mar(1));
end
M0=ma';
%********** find lambda **********
for k=1:tr
    latr=sprintf('[T%d]',k);
    xtr=findstr(s,latr);
    dels=findstr(s(xtr:length(s)),'TEMPORISATION');
    lambda(k)=getnum(s,xtr+dels(1));
    if lambda(k)==0 lambda(k)=1; end % Default value 1
end
lambda=lambda';
%********** find interconnections **********
xar=findstr(s,'[A');
p2t=findstr(s,'PLACE2TRANS');
sou=findstr(s,'SOURCE');
des=findstr(s,'DEST');
val=findstr(s,'VALUE');
Pre=zeros(pl,tr);
Post=zeros(pl,tr);
if ((size(xar,2)~=ar)||(size(p2t,2)~=ar)||(size(sou,2)~=ar)||(size(des,2)~=ar)||(size(val,2)~=ar))
    sprintf('Arcs red incorrectly')
end
for k=1:ar
    p2(k)=getnum(s,p2t(k));
    so(k)=getnum(s,sou(k))+1;
    de(k)=getnum(s,des(k))+1;
    va(k)=getnum(s,val(k));
    if p2(k)==1
        Pre(so(k),de(k))=va(k);
    else
        Post(de(k),so(k))=va(k);
    end
end
fclose(fid);


function [n] = getnum(str,i)
%gets a number from a string containing '='

%skip over the text
while strcmp(str(i),'=')==0
    i=i+1;
end
i=i+1;
j=i;	%number start position
while abs(str(i))~=13
    i=i+1;
end
n=str2num(str(j:i));
if isempty(n)
    n=0;
end
