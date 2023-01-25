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

function [p_semi,t_semi] = SimHPN_ptsemis(C)

p_semi=cal_semi(C)';  %Vectores columna
t_semi=cal_semi(C')';


function semi=cal_semi(C)

% Inicialización

cero=1e-3;
A = C;
D=eye(size(C,1));
colsC=size(C,2);

%D|A=[D,A]

for i=1:colsC,
    %
    % Buscar parejas de filas de D|A cuya combinación lineal positiva
    % anule la i-ésima columna de A
    %
    nA=zeros(1,size(A,2)); % Variables auxiliares para nuevas filas
    nD=zeros(1,size(D,2));
    rowsDA=size(A,1);
    for j=1:rowsDA-1,
        for k=j+1:rowsDA,
            if abs(A(k,i))>0
                sig=A(j,i)*A(k,i);
            else
                sig=1;
            end
            if sig < 0  % Anula la i-ésima columna con combinación positiva
                % nD y nA son las nuevas filas que hay que añadir a D|A
                nfD=D(j,:)*abs(A(k,i))+D(k,:)*abs(A(j,i));
                nD=[nD;nfD];
                nfA=A(j,:)*abs(A(k,i))+A(k,:)*abs(A(j,i));
                nA=[nA;nfA];
            end
        end
    end
    % Añadir las filas calculadas a D|A
    if size(nD,1)>1
        D=[D;nD(2:size(nD,1),:)];
        A=[A;nA(2:size(nA,1),:)];
    end

    %
    % Eliminar de D|A las filas en que la i-ésima columna de A sea no nula
    %
    rowsDA=size(A,1);
    nA=zeros(1,size(A,2)); % Variables auxiliares para nuevas filas
    nD=zeros(1,size(D,2));
    for j=1:rowsDA,
        if abs(A(j,i))<=cero    % Mantener fila j
            nfD=D(j,:);
            nD=[nD;nfD];
            nfA=A(j,:);
            nA=[nA;nfA];
        end
    end
    % Nuevos valores de D y A
    if size(nD,1)>1
        D=nD(2:size(nD,1),:);
        A=nA(2:size(nA,1),:);
    else
        D=zeros(1,size(nD,2));
        A=zeros(1,size(nA,2));
    end
end

% En este punto D es una matriz de enteros que contiene los semiflujos
% Pueden no ser mínimos y pueden estar repetidos

% Eliminación de semiflujos que contienen a otros

esse=ones(size(D,1),1); % 1 si la fila es semiflujo mínimo
supp=or(D,D);           % support de los semiflujos

for i=1:size(D,1)-1
    for j=i+1:size(D,1)
        uni=or(supp(i,:),supp(j,:));
        if min(supp(i,:)-uni)>=0  % i contiene a j
            esse(i,1)=0;
        elseif min(supp(j,:)-uni)>=0  % j contiene a i
            esse(j,1)=0; % En caso de ser iguales sólo se elimina uno
        end
    end
end

j=1;
for i=1:size(D,1)
    if esse(i,1)==1
        csemi(j,:)=D(i,:);
        j=j+1;
    end
end

% División de todos los semiflujos por el máximo común divisor

for i=1:size(csemi,1)
    mcd=csemi(i,1);
    for j=2:size(csemi,2)
        mcd=gcd(mcd,csemi(i,j));
    end
    if mcd>0
        semi(i,:)=csemi(i,:)/mcd;
    else
        semi(i,:)=csemi(i,:);
    end
end
