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


function [x,v,opt] = quadprog2(varargin)
% QUADPROG2 - Convex Quadratic Programming Solver
%             Featuring the SOLVOPT freeware optimizer
% 
%             New for version 1.1: * Significant Speed Improvement
%                                  * Geometric Preconditioning
%                                  * Improved Error Checking
% 
% USAGE:
% [x,v]     = quadprog2(H,f,A,b)
% [x,v]     = quadprog2(H,f,A,b,guess)
% [x,v,opt] = ...
% 
% Minimizes the function v = 0.5*x'*H*x + f*x
% subject to the constraint A*x <= b.
% Initial guess is optional.
% 
% ("opt" returns SOLVOPT data for advanced use. Details are available in
% the SOLVOPT documentation at the website identified below.)
% 
% Notes: (1) For a problem with 100 variables and 300 constraints, you will
%            often get a result in under 5 seconds. However, sometimes
%            the optimizer has to work longer (see below) for difficult
%            optimizations. Alerts are provided. (Note: The calculation
%            time is more sensitive to the number of variables than it is
%            to the number of contraints.)
%        (2) Geometric preconditioning is undertaken for 10 or more
%            dimensions to greatly reduce calculation time. (With fewer
%            than 10 dimensions, there is negligible benefit, so the
%            preconditioning calculations are omitted.)
%        (3) Geometric preconditioning can impair the convergence of some
%            difficult optimizations. When this occurs, the optimization
%            is attempted again without the preconditioning.
%        (4) x and guess are column vectors. f is a row vector.
%            They will be converted if necessary.
%        (5) This m-file incorporates the SOLVOPT version 1.1 freeware
%            optimizer, which has been wholly reproduced, except for
%            a few slight modifications for convenience in parameter passing.
%        (6) SolvOpt is a general nonlinear local optimizer,
%            written by Alexei Kuntsevich & Franz Kappel, and
%            is available (as of this writing) as freeware from:
%            http://www.uni-graz.at/imawww/kuntsevich/solvopt/
%        (7) This Matlab function requires a convex QP problem
%            with a positive-definite symmetric matrix H.
%            This is a somewhat trivial application of
%            a general solver like SOLVOPT, but the use of precomputed
%            gradient vectors herein makes the solution fast enough
%            to warrant use.
%        (8) Any local solution of a convex QP is also a global solution.
%            Hence, your results will be globally optimal.
%        (9) Relative precision in the objective function is set to 1e-6.
%       (10) Absolute precision in constraint violation is 1e-6 or better.
%       (11) This program does not require the Optimization Toolbox
%       (12) ver 1.0: intial writing, Michael Kleder, June 2005
%       (13) ver 1.1: geometric preconditioning, Michael Kleder, July 2005
% 
% EXAMPLE:
% % Convex QP with 100 variables and 300 constraints:
% n = 100;
% c = 300;
% H = rand(n);
% H=H*H';
% f=rand(1,n);
% A=rand(c,n)*2-1;
% b=ones(c,1);
% tic
% x = quadprog2(H,f,A,b)
% toc

[H,f,A,b]=deal(varargin{1:4});
if any(eig(H)<=0)
    error('Non-convex. (Quadratic form is not positive definite.)')
end
if any(H~=H')
    error('Quadratic form must be symmetric.')
end
f=f(:)';
n=size(H,2);
% try a quick check on unconstrained solution
HI=inv(H);
mqf=-HI*f'; % inverse exists since H is positive definite
if all(A*mqf<=b(:))
    x=mqf;
    v = .5*x'*H*x + f*x;
    opt=[];
    return
end
if nargin < 5
    x=mqf+(rand(n,1)-.5);
else
    x=varargin{5};
    x=x(:);
end

if n>=10
    % geometric pre-conditioning:
    bn = b - A*mqf;
    fn = 0*f;
    T = sqrtm(HI);
    An = A*T;
    nA=sqrt(sum(An.^2,2));
    An = An ./ repmat(nA,[1 size(An,2)]);
    bn = bn ./ nA;
    Hn = eye(size(H));
    % invoke nonlinear solver:
    pass={Hn,fn,An,bn};
    %opt=[-1,1e-4,1e-6,15000,0,1e-8,2.5,1e-11]; % defaults
    opt=[-1,1e-6,1e-8,15e5,0,1e-8,2.5,1e-11]; % increase precision
    %opt=[-1,1e-14,1e-14,15e5,0,1e-14,2.5,1e-12]; % very high precision
    %[x,v,opt] = solvopt2(x,@obj,[],opt,[],[],pass);
    opt(6) = opt(6) * .1; % precondition needs constraint tolerance reduced
    [x,v,opt] = solvopt2(x,@obj,@gobj,opt,@pen,@gpen,pass);
    opt(6) = opt(6) * .1;
    x = T*x+mqf;
    v = .5*x'*H*x + f*x;
    if max(A*x-b) <= 1e-6
        return
    else
        disp('Difficult optimization. Omitting geometric preconditioning.')
        pass={H,f,A,b};
        [x,v,opt] = solvopt2(x,@obj,@gobj,opt,@pen,@gpen,pass);
        if max(A*x-b) > opt(6)
            disp('Feasible optimum not found. Searching again.')
            x=mqf+(rand(n,1)-.5);
            [x,v,opt] = solvopt2(x,@obj,@gobj,opt,@pen,@gpen,pass);
            if max(A*x-b) > opt(6)
                error('Optimization failure: Feasible optimum not found.')
            end
        end
        return
    end
else % low-dimensional
    pass={H,f,A,b};
    %opt=[-1,1e-4,1e-6,15000,0,1e-8,2.5,1e-11]; % defaults
    opt=[-1,1e-6,1e-8,15e5,0,1e-8,2.5,1e-11]; % increase precision
    %opt=[-1,1e-14,1e-14,15e5,0,1e-14,2.5,1e-12]; % very high precision
    [x,v,opt] = solvopt2(x,@obj,@gobj,opt,@pen,@gpen,pass);
    if max(A*x-b) > opt(6)
        disp('Feasible optimum not found. Searching again.')
        x=mqf+(rand(n,1)-.5);
        [x,v,opt] = solvopt2(x,@obj,@gobj,opt,@pen,@gpen,pass);
        if max(A*x-b) > opt(6)
            error('Optimization failure: feasible optimum not found.')
        end
    end
end
return
function q = obj(x,pass) % objective function
H=pass{1};
f=pass{2};
q = .5*x'*H*x+f*x;
return
function q = gobj(x,pass); % gradient of objective function
H=pass{1};
f=pass{2};
q = H*x+f';
return
function q=pen(x,pass) % contraint penalty
A=pass{3};
b=pass{4};
q=max(0,A*x-b);
q=sum(q(:));
return
function q=gpen(x,pass) % gradient of contraint penalty
A=pass{3};
b=pass{4};
q=A*x-b;
q(q<0)=0;
q = q'*A/sqrt(q'*q);
return
function [x,f,options]=solvopt2(x,fun,grad,options,func,gradc,pass)
% added 'pass' input to pass params to fun/func/grad/gradc; also quieted warnings
% Usage:
% [x,f,options]=solvopt(x,fun,grad,options,func,gradc)
% The function SOLVOPT performs a modified version of Shor's r-algorithm in
% order to find a local minimum resp. maximum of a nonlinear function
% defined on the n-dimensional Euclidean space
% or
% a solution of a nonlinear constrained problem:
% min { f(x): g(x) (<)= 0, g(x) in R(m), x in R(n) }
% Arguments:
% x       is the n-vector (row or column) of the coordinates of the starting
%         point,
% fun     is the name of an M-file (M-function) which computes the value
%         of the objective function <fun> at a point x,
%         synopsis: f=fun(x)
% grad    is the name of an M-file (M-function) which computes the gradient
%         vector of the function <fun> at a point x,
%         synopsis: g=grad(x)
% func    is the name of an M-file (M-function) which computes the MAXIMAL
%         RESIDUAL(!) for a set of constraints at a point x,
%         synopsis: fc=func(x)
% gradc   is the name of an M-file (M-function) which computes the gradient
%         vector for the maximal residual consyraint at a point x,
%         synopsis: gc=gradc(x)
% options is a row vector of optional parameters:
%    options(1)= H, where sign(H)=-1 resp. sign(H)=+1 means minimize
%        resp. maximize <fun> (valid only for unconstrained problem)
%        and H itself is a factor for the initial trial step size
%        (options(1)=-1 by default),
%    options(2)= relative error for the argument in terms of the
%        infinity-norm (1.e-4 by default),
%    options(3)= relative error for the function value (1.e-6 by default),
%    options(4)= limit for the number of iterations (15000 by default),
%    options(5)= control of the display of intermediate results and error
%        resp. warning messages (default value is 0, i.e., no intermediate
%        output but error and warning messages, see more in the manual),
%    options(6)= admissible maximal residual for a set of constraints
%        (options(6)=1e-8 by default, see more in the manual),
%   *options(7)= the coefficient of space dilation (2.5 by default),
%   *options(8)= lower bound for the stepsize used for the difference
%        approximation of gradients (1e-12 by default, see more in the manual).
%  (* ... changes should be done with care)
% Returned values:
% x       is the optimizer (row resp. column),
% f       is the optimum function value,
% options returns the values of the counters
%    options(9),  the number of iterations, if positive,
%        or an abnormal stop code, if negative (see more in the manual),
%    options(10), the number of objective
%    options(11), the number of gradient evaluations,
%    options(12), the number of constraint function evaluations,
%    options(13), the number of constraint gradient evaluations.
% ____________________________________________________________________________

shhh = 1; % silence warnings (not errors)

% strings: ----{
errmes='SolvOpt error:';
wrnmes='SolvOpt warning:';
error1='No function name and/or starting point passed to the function.';
error2='Argument X has to be a row or column vector of dimension > 1.';
error30='<fun> returns an empty string.';
error31='Function value does not exist (NaN is returned).';
error32='Function equals infinity at the point.';
error40='<grad> returns an improper matrix. Check the dimension.';
error41='Gradient does not exist (NaN is returned by <grad>).';
error42='Gradient equals infinity at the starting point.';
error43='Gradient equals zero at the starting point.';
error50='<func> returns an empty string.';
error51='<func> returns NaN at the point.';
error52='<func> returns infinite value at the point.';
error60='<gradc> returns an improper vector. Check the dimension';
error61='<gradc> returns NaN at the point.';
error62='<gradc> returns infinite vector at the point.';
error63='<gradc> returns zero vector at an infeasible point.';
error5='Function is unbounded.';
error6='Choose another starting point.';
warn1= 'Gradient is zero at the point, but stopping criteria are not fulfilled.';
warn20='Normal re-setting of a transformation matrix.' ;
warn21='Re-setting due to the use of a new penalty coefficient.' ;
warn4= 'Iterations limit exceeded.';
warn31='The function is flat in certain directions.';
warn32='Trying to recover by shifting insensitive variables.';
warn09='Re-run from recorded point.';
warn08='Ravine with a flat bottom is detected.';
termwarn0='SolvOpt: Normal termination.';
termwarn1='SolvOpt: Termination warning:';
appwarn='The above warning may be reasoned by inaccurate gradient approximation';
endwarn=[...
    'Premature stop is possible. Try to re-run the routine from the obtained point.               ';...
    'Result may not provide the optimum. The function apparently has many extremum points.        ';...
    'Result may be inaccurate in the coordinates. The function is flat at the optimum.            ';...
    'Result may be inaccurate in a function value. The function is extremely steep at the optimum.'];
% ----}

% ARGUMENTS PASSED ----{
if nargin<2           % Function and/or starting point are not specified
    options(9)=-1; disp(errmes);  disp(error1);   return
end
if nargin<3,   app=1;             % No user-supplied gradients
elseif isempty(grad),  app=1;
else,  app=0;                     % Exact gradients are supplied
end

% OPTIONS ----{
doptions=[-1,1.e-4,1.e-6,15000,0,1.e-8,2.5,1e-11];
if nargin<4,  options=doptions;
elseif isempty(options), options=doptions;
else,
    % Replace default options by user specified options:
    ii=find(options~=0);doptions(ii)=options(ii);
    options=doptions;
end
% Check the values:
options([2:4,6:8])=abs(options([2:4,6:8]));
options(2:3)=max(options(2:3),[1.e-12,1.e-12]);
options(2)=max(options(8)*1.e2,options(2));
options(2:3)=min(options(2:3),[1,1]);
options(6)=max(options(6),1e-12);
options(7)=max([options(7),1.5]);
options(8)=max(options(8),1e-11);
% ----}

if nargin<5,   constr=0;          % Unconstrained problem
elseif isempty(func),  constr=0;
else,  constr=1;                  % Constrained problem
    if nargin<6, appconstr=1; t=3; % No user-supplied gradients for constraints
    elseif isempty(gradc),
        appconstr=1;
    else,  appconstr=0;            % Exact gradients of constraints are supplied
    end
end
% ----}

% STARTING POINT ----{
if max(size(x))<=1,      disp(errmes);  disp(error2);
    options(9)=-2; return
elseif size(x,2)==1,     n=size(x,1);  x=x'; trx=1;
elseif size(x,1)==1,     n=size(x,2);        trx=0;
else,                    disp(errmes);  disp(error2);
    options(9)=-2; return
end
% ----}

% WORKING CONSTANTS AND COUNTERS ----{

options(10)=0; options(11)=0;      % function and gradient calculations
if constr
    options(12)=0; options(13)=0;      % same for constraints
end
epsnorm=1.e-15;epsnorm2=1.e-30;    % epsilon & epsilon^2

if constr, h1=-1;                  % NLP: restricted to minimization
    cnteps=options(6);                % Max. admissible residual
else, h1=sign(options(1));         % Minimize resp. maximize a function
end

k=0;                               % Iteration counter

wdef=1/options(7)-1;               % Default space transf. coeff.

%Gamma control ---{
ajb=1+.1/n^2;                    % Base I
ajp=20;
ajpp=ajp;                        % Start value for the power
ajs=1.15;                        % Base II
knorms=0; gnorms=zeros(1,10);    % Gradient norms stored
%---}

%Display control ---{
if options(5)<=0, dispdata=0;
    if options(5)==-1, dispwarn=0; else, dispwarn=1; end
else, dispdata=round(options(5)); dispwarn=1;
end,  ld=dispdata;
%---}

%Stepsize control ---{
dq=5.1;                          % Step divider (at f_{i+1}>gamma*f_{i})
du20=2;du10=1.5;du03=1.05;       % Step multipliers (at certain steps made)
kstore=3;nsteps=zeros(1,kstore); % Steps made at the last 'kstore' iterations
if app, des=6.3;                 % Desired number of steps per 1-D search
else,   des=3.3; end
mxtc=3;                          % Number of trial cycles (steep wall detect)
%---}
termx=0; limxterm=50;              % Counter and limit for x-criterion

ddx   =max(1e-11,options(8));      % stepsize for gradient approximation

low_bound=-1+1e-4;                 % Lower bound cosine used to detect a ravine

ZeroGrad=n*1.e-16;                 % Lower bound for a gradient norm

nzero=0;                           % Zero-gradient events counter
% Lower bound for values of variables taking into account
lowxbound=max([options(2),1e-3]);
% Lower bound for function values to be considered as making difference
lowfbound=options(3)^2;
krerun=0;                          % Re-run events counter
detfr=options(3)*100;              % relative error for f/f_{record}
detxr=options(2)*10;               % relative error for norm(x)/norm(x_{record})

warnno=0;                          % the number of warn.mess. to end with

kflat=0;                           % counter for points of flatness
stepvanish=0;                      % counter for vanished steps
stopf=0;
% ----}  End of setting constants
% ----}  End of the preamble

% COMPUTE THE FUNCTION  ( FIRST TIME ) ----{
if trx,  f=feval(fun,x',pass);
else,    f=feval(fun,x,pass);  end
options(10)=options(10)+1;
if isempty(f),      if dispwarn,disp(errmes);disp(error30);end
    options(9)=-3; if trx, x=x';end, return
elseif isnan(f),    if dispwarn,disp(errmes);disp(error31);disp(error6);end
    options(9)=-3; if trx, x=x';end, return
elseif abs(f)==Inf, if dispwarn,disp(errmes);disp(error32);disp(error6);end
    options(9)=-3; if trx, x=x';end, return
end
xrec=x; frec=f;     % record point and function value
% Constrained problem
if constr,  fp=f; kless=0;
    if trx,  fc=feval(func,x',pass);
    else,    fc=feval(func,x,pass);  end
    if isempty(fc),
        if dispwarn,disp(errmes);disp(error50);end
        options(9)=-5; if trx, x=x';end, return
    elseif isnan(fc),
        if dispwarn,disp(errmes);disp(error51);disp(error6);end
        options(9)=-5; if trx, x=x';end, return
    elseif abs(fc)==Inf,
        if dispwarn,disp(errmes);disp(error52);disp(error6);end
        options(9)=-5; if trx, x=x';end, return
    end
    options(12)=options(12)+1;
    PenCoef=1;                              % first rough approximation
    if fc<=cnteps,  FP=1; fc=0;             % feasible point
    else,           FP=0;                   % infeasible point
    end
    f=f+PenCoef*fc;
end
% ----}
% COMPUTE THE GRADIENT ( FIRST TIME ) ----{
if app,   deltax=h1*ddx*ones(size(x));
    if constr, if trx, g=apprgrdn2(x',fp,fun,pass,deltax',1);
        else,   g=apprgrdn2(x ,fp,fun,pass,deltax,1); end
    else,      if trx, g=apprgrdn2(x',f,fun,pass,deltax',1);
        else,   g=apprgrdn2(x ,f,fun,pass,deltax,1); end
    end, options(10)=options(10)+n;
else,     if trx,  g=feval(grad,x',pass);
    else,    g=feval(grad,x,pass);   end
    options(11)=options(11)+1;
end
if size(g,2)==1, g=g'; end, ng=norm(g);
if size(g,2)~=n,    if dispwarn,disp(errmes);disp(error40);end
    options(9)=-4; if trx, x=x';end, return
elseif isnan(ng),   if dispwarn,disp(errmes);disp(error41);disp(error6);end
    options(9)=-4; if trx, x=x';end, return
elseif ng==Inf,     if dispwarn,disp(errmes);disp(error42);disp(error6);end
    options(9)=-4; if trx, x=x';end, return
elseif ng<ZeroGrad, if dispwarn,disp(errmes);disp(error43);disp(error6);end
    options(9)=-4; if trx, x=x';end, return
end
if constr, if ~FP
        if appconstr,
            deltax=sign(x); idx=find(deltax==0);
            deltax(idx)=ones(size(idx));  deltax=ddx*deltax;
            if trx, gc=apprgrdn2(x',fc,func,pass,deltax',0);
            else,   gc=apprgrdn2(x ,fc,func,pass,deltax ,0); end
            options(12)=options(12)+n;
        else,     if trx,  gc=feval(gradc,x',pass);
            else,    gc=feval(gradc,x,pass); end
            options(13)=options(13)+1;
        end
        if size(gc,2)==1, gc=gc'; end, ngc=norm(gc);
        if size(gc,2)~=n,
            if dispwarn,disp(errmes);disp(error60);end
            options(9)=-6; if trx, x=x';end, return
        elseif isnan(ngc),
            if dispwarn,disp(errmes);disp(error61);disp(error6);end
            options(9)=-6; if trx, x=x';end, return
        elseif ngc==Inf,
            if dispwarn,disp(errmes);disp(error62);disp(error6);end
            options(9)=-6; if trx, x=x';end, return
        elseif ngc<ZeroGrad,
            if dispwarn,disp(errmes);disp(error63);end
            options(9)=-6; if trx, x=x';end, return
        end
        g=g+PenCoef*gc; ng=norm(g);
    end, end
grec=g; nng=ng;
% ----}
% INITIAL STEPSIZE
h=h1*sqrt(options(2))*max(abs(x));     % smallest possible stepsize
if abs(options(1))~=1,
    h=h1*max(abs([options(1),h]));     % user-supplied stepsize
else,
    h=h1*max(1/log(ng+1.1),abs(h));    % calculated stepsize
end

% RESETTING LOOP ----{
while 1,
    kcheck=0;                        % Set checkpoint counter.
    kg=0;                            % stepsizes stored
    kj=0;                            % ravine jump counter
    B=eye(n);                        % re-set transf. matrix to identity
    fst=f; g1=g;  dx=0;
    % ----}

    % MAIN ITERATIONS ----{

    while 1,
        k=k+1;kcheck=kcheck+1;
        laststep=dx;

        % ADJUST GAMMA --{
        gamma=1+max([ajb^((ajp-kcheck)*n),2*options(3)]);
        gamma=min([gamma,ajs^max([1,log10(nng+1)])]);
        % --}
        gt=g*B;   w=wdef;
        % JUMPING OVER A RAVINE ----{
        if (gt/norm(gt))*(g1'/norm(g1))<low_bound
            if kj==2, xx=x;  end,  if kj==0, kd=4;  end,
            kj=kj+1;  w=-.9; h=h*2;             % use large coef. of space dilation
            if kj>2*kd,     kd=kd+1;  warnno=1;
                if any(abs(x-xx)<epsnorm*abs(x)), % flat bottom is detected
                    if dispwarn & ~shhh,disp(wrnmes);disp(warn08); end
                end
            end
        else, kj=0;
        end
        % ----}
        % DILATION ----{
        z=gt-g1;
        nrmz=norm(z);
        if(nrmz>epsnorm*norm(gt))
            z=z/nrmz;
            g1=gt+w*(z*gt')*z;  B=B+w*(B*z')*z;
        else
            z=zeros(1,n); nrmz=0; g1=gt;
        end
        d1=norm(g1);  g0=(g1/d1)*B';
        % ----}
        % RESETTING ----{
        if kcheck>1
            idx=find(abs(g)>ZeroGrad); numelem=size(idx,2);
            if numelem>0, grbnd=epsnorm*numelem^2;
                if all(abs(g1(idx))<=abs(g(idx))*grbnd) | nrmz==0
                    if dispwarn & ~shhh,  disp(wrnmes);  disp(warn20); end
                    if abs(fst-f)<abs(f)*.01, ajp=ajp-10*n;
                    else, ajp=ajpp; end
                    h=h1*dx/3; k=k-1;
                    break
                end
            end
        end
        % ----}
        % STORE THE CURRENT VALUES AND SET THE COUNTERS FOR 1-D SEARCH
        xopt=x;fopt=f;   k1=0;k2=0;ksm=0;kc=0;knan=0;  hp=h;
        if constr, Reset=0; end
        % 1-D SEARCH ----{
        while 1,
            x1=x;f1=f;
            if constr, FP1=FP; fp1=fp; end
            x=x+hp*g0;
            % FUNCTION VALUE
            if trx, f=feval(fun,x',pass);
            else,   f=feval(fun,x,pass );  end
            options(10)=options(10)+1;
            if h1*f==Inf
                if dispwarn, disp(errmes); disp(error5); end
                options(9)=-7; if trx, x=x';end, return
            end
            if constr, fp=f;
                if trx,fc=feval(func,x',pass);
                else,  fc=feval(func,x,pass);end
                options(12)=options(12)+1;
                if  isnan(fc),
                    if dispwarn,disp(errmes);disp(error51);disp(error6);end
                    options(9)=-5; if trx, x=x';end, return
                elseif abs(fc)==Inf,
                    if dispwarn,disp(errmes);disp(error52);disp(error6);end
                    options(9)=-5; if trx, x=x';end, return
                end
                if fc<=cnteps,   FP=1; fc=0;
                else,            FP=0;
                    fp_rate=(fp-fp1);
                    if fp_rate<-epsnorm
                        if ~FP1
                            PenCoefNew=-15*fp_rate/norm(x-x1);
                            if PenCoefNew>1.2*PenCoef,
                                PenCoef=PenCoefNew; Reset=1; kless=0; f=f+PenCoef*fc; break
                            end
                        end
                    end
                end
                f=f+PenCoef*fc;
            end
            if abs(f)==Inf | isnan(f)
                if dispwarn & ~shhh, disp(wrnmes);
                    if isnan(f), disp(error31); else, disp(error32); end
                end
                if ksm | kc>=mxtc, options(9)=-3; if trx, x=x';end, return
                else, k2=k2+1;k1=0; hp=hp/dq; x=x1;f=f1; knan=1;
                    if constr, FP=FP1; fp=fp1; end
                end
                % STEP SIZE IS ZERO TO THE EXTENT OF EPSNORM
            elseif all(abs(x-x1)<abs(x)*epsnorm),
                stepvanish=stepvanish+1;
                if stepvanish>=5,
                    options(9)=-14;
                    if dispwarn, disp(termwarn1);
                        disp(endwarn(4,:)); end
                    if trx,x=x';end,  return
                else, x=x1; f=f1; hp=hp*10; ksm=1;
                    if constr, FP=FP1; fp=fp1; end
                end
                % USE SMALLER STEP
            elseif h1*f<h1*gamma^sign(f1)*f1
                if ksm,break,end,  k2=k2+1;k1=0; hp=hp/dq; x=x1;f=f1;
                if constr, FP=FP1; fp=fp1; end
                if kc>=mxtc, break, end
                % 1-D OPTIMIZER IS LEFT BEHIND
            else   if h1*f<=h1*f1, break,  end
                % USE LARGER STEP
                k1=k1+1; if k2>0, kc=kc+1; end, k2=0;
                if k1>=20,      hp=du20*hp;
                elseif k1>=10,  hp=du10*hp;
                elseif k1>=3,   hp=du03*hp;
                end
            end
        end
        % ----}  End of 1-D search
        % ADJUST THE TRIAL STEP SIZE ----{
        dx=norm(xopt-x);
        if kg<kstore,  kg=kg+1;  end
        if kg>=2,  nsteps(2:kg)=nsteps(1:kg-1); end
        nsteps(1)=dx/(abs(h)*norm(g0));
        kk=sum(nsteps(1:kg).*[kg:-1:1])/sum([kg:-1:1]);
        if     kk>des, if kg==1,  h=h*(kk-des+1);
            else,   h=h*sqrt(kk-des+1); end
        elseif kk<des,         h=h*sqrt(kk/des);
        end

        stepvanish=stepvanish+ksm;
        % ----}
        % COMPUTE THE GRADIENT ----{
        if app,
            deltax=sign(g0); idx=find(deltax==0);
            deltax(idx)=ones(size(idx));  deltax=h1*ddx*deltax;
            if constr,  if trx,  g=apprgrdn2(x',fp,fun,pass,deltax',1);
                else,    g=apprgrdn2(x ,fp,fun,pass,deltax ,1);    end
            else,       if trx,  g=apprgrdn2(x',f,fun,pass,deltax',1);
                else,    g=apprgrdn2(x ,f,fun,pass,deltax ,1);    end
            end,  options(10)=options(10)+n;
        else
            if trx,  g=feval(grad,x',pass);
            else,    g=feval(grad,x ,pass); end
            options(11)=options(11)+1;
        end
        if size(g,2)==1, g=g'; end,    ng=norm(g);
        if isnan(ng),
            if dispwarn, disp(errmes); disp(error41); end
            options(9)=-4; if trx, x=x'; end, return
        elseif ng==Inf,
            if dispwarn,disp(errmes);disp(error42);end
            options(9)=-4; if trx, x=x';end, return
        elseif ng<ZeroGrad,
            if dispwarn & ~shhh,disp(wrnmes);disp(warn1);end
            ng=ZeroGrad;
        end
        % Constraints:
        if constr, if ~FP
                if ng<.01*PenCoef
                    kless=kless+1;
                    if kless>=20, PenCoef=PenCoef/10; Reset=1; kless=0; end
                else, kless=0;
                end
                if appconstr,
                    deltax=sign(x); idx=find(deltax==0);
                    deltax(idx)=ones(size(idx));  deltax=ddx*deltax;
                    if trx, gc=apprgrdn2(x',fc,func,pass,deltax',0);
                    else,   gc=apprgrdn2(x ,fc,func,pass,deltax ,0); end
                    options(12)=options(12)+n;
                else, if trx,  gc=feval(gradc,x',pass);
                    else,    gc=feval(gradc,x ,pass); end
                    options(13)=options(13)+1;
                end
                if size(gc,2)==1, gc=gc'; end, ngc=norm(gc);
                if     isnan(ngc),
                    if dispwarn,disp(errmes);disp(error61);end
                    options(9)=-6; if trx, x=x';end, return
                elseif ngc==Inf,
                    if dispwarn,disp(errmes);disp(error62);end
                    options(9)=-6; if trx, x=x';end, return
                elseif ngc<ZeroGrad & ~appconstr,
                    if dispwarn,disp(errmes);disp(error63);end
                    options(9)=-6; if trx, x=x';end, return
                end
                g=g+PenCoef*gc; ng=norm(g);
                if Reset, if dispwarn & ~shhh,  disp(wrnmes);  disp(warn21); end
                    h=h1*dx/3; k=k-1; nng=ng; break
                end
            end, end
        if h1*f>h1*frec, frec=f; xrec=x; grec=g; end
        % ----}
        if ng>ZeroGrad,
            if knorms<10,  knorms=knorms+1;  end
            if knorms>=2,  gnorms(2:knorms)=gnorms(1:knorms-1); end
            gnorms(1)=ng;
            nng=(prod(gnorms(1:knorms)))^(1/knorms);
        end

        % DISPLAY THE CURRENT VALUES ----{
        if k==ld
            disp('Iter.# ..... Function ... Step Value ... Gradient Norm ');
            disp(sprintf('%5i   %13.5e   %13.5e     %13.5e',k,f,dx,ng));
            ld=k+dispdata;
        end
        %----}
        % CHECK THE STOPPING CRITERIA ----{
        termflag=1;
        if constr, if ~FP, termflag=0; end, end
        if kcheck<=5, termflag=0; end
        if knan, termflag=0; end
        if kc>=mxtc, termflag=0; end
        % ARGUMENT
        if termflag
            idx=find(abs(x)>=lowxbound);
            if isempty(idx) | all(abs(xopt(idx)-x(idx))<=options(2)*abs(x(idx)))
                termx=termx+1;
                % FUNCTION
                if abs(f-frec)> detfr * abs(f)    & ...
                        abs(f-fopt)<=options(3)*abs(f) & ...
                        krerun<=3                      & ...
                        ~constr
                    if any(abs(xrec(idx)-x(idx))> detxr * abs(x(idx)))
                        if dispwarn & ~shhh,disp(wrnmes);disp(warn09);end
                        x=xrec; f=frec; g=grec; ng=norm(g); krerun=krerun+1;
                        h=h1*max([dx,detxr*norm(x)])/krerun;
                        warnno=2; break
                    else, h=h*10;
                    end
                elseif  abs(f-frec)> options(3)*abs(f)    & ...
                        norm(x-xrec)<options(2)*norm(x) & constr

                elseif  abs(f-fopt)<=options(3)*abs(f)  | ...
                        abs(f)<=lowfbound               | ...
                        (abs(f-fopt)<=options(3) & termx>=limxterm )
                    if stopf
                        if dx<=laststep
                            if warnno==1 & ng<sqrt(options(3)), warnno=0; end
                            if ~app, if any(abs(g)<=epsnorm2), warnno=3; end, end
                            if warnno~=0, options(9)=-warnno-10;
                                if dispwarn, disp(termwarn1);
                                    disp(endwarn(warnno,:));
                                    if app, disp(appwarn); end
                                end
                            else, options(9)=k;
                                if dispwarn, disp(termwarn0); end
                            end
                            if trx,x=x';end,  return
                        end
                    else, stopf=1;
                    end
                elseif dx<1.e-12*max(norm(x),1) & termx>=limxterm
                    options(9)=-14;
                    if dispwarn, disp(termwarn1); disp(endwarn(4,:));
                        if app, disp(appwarn); end
                    end
                    x=xrec; f=frec;
                    if trx,x=x';end,  return
                else, stopf=0;
                end
            end
        end
        % ITERATIONS LIMIT
        if(k==options(4))
            options(9)=-9; if trx, x=x'; end,
            if dispwarn & ~shhh, disp(wrnmes);  disp(warn4); end
            return
        end
        % ----}
        % ZERO GRADIENT ----{
        if constr
            if ng<=ZeroGrad,
                if dispwarn,  disp(termwarn1);  disp(warn1); end
                options(9)=-8; if trx,x=x';end,return
            end
        else
            if ng<=ZeroGrad,        nzero=nzero+1;
                if dispwarn & ~shhh, disp(wrnmes);  disp(warn1);  end
                if nzero>=3,  options(9)=-8; if trx,x=x';end,return, end
                g0=-h*g0/2;
                for i=1:10,
                    x=x+g0;
                    if trx, f=feval(fun,x',pass);
                    else,   f=feval(fun,x ,pass); end
                    options(10)=options(10)+1;
                    if abs(f)==Inf
                        if dispwarn, disp(errmes);  disp(error32);  end
                        options(9)=-3;if trx,x=x';end,return
                    elseif isnan(f),
                        if dispwarn, disp(errmes);  disp(error32);  end
                        options(9)=-3;if trx,x=x';end,return
                    end
                    if app,
                        deltax=sign(g0); idx=find(deltax==0);
                        deltax(idx)=ones(size(idx));  deltax=h1*ddx*deltax;
                        if trx,  g=apprgrdn2(x',f,fun,pass,deltax',1);
                        else,    g=apprgrdn2(x ,f,fun,pass,deltax ,1);    end
                        options(10)=options(10)+n;
                    else
                        if trx,  g=feval(grad,x',pass);
                        else,    g=feval(grad,x ,pass);   end
                        options(11)=options(11)+1;
                    end
                    if size(g,2)==1, g=g'; end,       ng=norm(g);
                    if ng==Inf
                        if dispwarn, disp(errmes);  disp(error42); end
                        options(9)=-4; if trx, x=x'; end, return
                    elseif isnan(ng)
                        if dispwarn, disp(errmes);  disp(error41); end
                        options(9)=-4; if trx, x=x'; end, return
                    end
                    if ng>ZeroGrad, break, end
                end
                if ng<=ZeroGrad,
                    if dispwarn,  disp(termwarn1);  disp(warn1); end
                    options(9)=-8; if trx,x=x';end,return
                end
                h=h1*dx; break
            end
        end
        % ----}
        % FUNCTION IS FLAT AT THE POINT ----{
        if ~constr & abs(f-fopt)<abs(fopt)*options(3) & kcheck>5 & ng<1
            idx=find(abs(g)<=epsnorm2); ni=size(idx,2);
            if ni>=1 & ni<=n/2 & kflat<=3, kflat=kflat+1;
                if dispwarn & ~shhh,  disp(wrnmes); disp(warn31); end, warnno=1;
                x1=x; fm=f;
                for j=idx, y=x(j); f2=fm;
                    if y==0, x1(j)=1; elseif abs(y)<1, x1(j)=sign(y); else, x1(j)=y; end
                    for i=1:20, x1(j)=x1(j)/1.15;
                        if trx, f1=feval(fun,x1',pass);
                        else,   f1=feval(fun,x1 ,pass); end
                        options(10)=options(10)+1;
                        if abs(f1)~=Inf & ~isnan(f1),
                            if h1*f1>h1*fm,     y=x1(j); fm=f1;
                            elseif h1*f2>h1*f1, break
                            elseif f2==f1,      x1(j)=x1(j)/1.5;
                            end, f2=f1;
                        end
                    end
                    x1(j)=y;
                end
                if h1*fm>h1*f
                    if app,
                        deltax=h1*ddx*ones(size(deltax));
                        if trx,  gt=apprgrdn2(x1',fm,fun,pass,deltax',1);
                        else,    gt=apprgrdn2(x1 ,fm,fun,pass,deltax ,1);    end
                        options(10)=options(10)+n;
                    else
                        if trx,  gt=feval(grad,x1',pass);
                        else,    gt=feval(grad,x1 ,pass); end
                        options(11)=options(11)+1;
                    end
                    if size(gt,2)==1, gt=gt'; end,       ngt=norm(gt);
                    if ~isnan(ngt) & ngt>epsnorm2,
                        if dispwarn,  disp(warn32); end
                        options(3)=options(3)/5;
                        x=x1; g=gt; ng=ngt; f=fm; h=h1*dx/3; break
                    end
                end
            end
        end
        % ----}
    end % iterations
end % restart
% end of the function
%

function g = apprgrdn2(x,f,fun,pass,deltax,obj)
% added 'pass' input to pass params to called function
% Usage:
% g = apprgrdn(x,f,fun,deltax,obj)
% Function apprgrdn.m performs the finite difference approximation
% of the gradient <g> at a point <x>.
% <f> is the calculated function value at a point <x>,
% <fun> is the name of the Matlab function, which calculates function values
% <deltax> is a vector of the relative stepsizes,
% <obj> is the flag indicating whether the gradient of the objective
%        function (1) or the constraint function (0) is to be calculated.
%
n=max(size(x)); ee=ones(size(x));
di=abs(x); idx=find(di<5e-15); di(idx)=5e-15*ee(idx);
di=deltax.*di;
if obj, idx=find(abs(di)<2e-10); di(idx)=2e-10*sign(di(idx));
else,   idx=find(abs(di)<5e-15); di(idx)=5e-15*sign(di(idx));
end
y=x;
for i=1:n
    y(i)=x(i)+di(i);
    fi=feval(fun,y,pass);
    if obj, if fi==f,
            for j=1:3
                di(i)=di(i)*10;  y(i)=x(i)+di(i);
                fi=feval(fun,y,pass); if fi~=f, break, end
            end
        end, end
    g(i)=(fi-f)/di(i);
    if obj, if any(idx==i)
            y(i)=x(i)-di(i);
            fi=feval(fun,y,pass);
            g(i)=.5*(g(i)+(f-fi)/di(i));
        end, end
    y(i)=x(i);
end

