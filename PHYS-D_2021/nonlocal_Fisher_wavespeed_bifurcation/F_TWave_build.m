function [] = F_TWave_build(ORDER)
% This function builds realified versions of function_F.m, its
% Frechet derivative and the quadratic Taylor expansion of this derivative
% for the radii polynomial. This is done using symbolics. These functions 
% files are output to the folder \Data. The ORDER input specifies to what
% order to Taylor expand DF. ORDER=2 will kill off every polynomial term.
if nargin==0
    ORDER=2;
end
addpath('Data')
X = sym('X',[2,1],'real');
Xprime = sym('Xprime',[2,1],'real');
cpseed = sym('cspeed','real');
omega = sym('omega','real');
v = sym('v',[2,1],'real');
lamprime = sym('lamprime',[2,1],'real');
vprime = sym('vprime',[2,1],'real');
para = sym('para',[3,1],'real');

[F1,F2,F3,F4,~] = function_F(X,Xprime,cpseed,omega,v(1)+1i*v(2),...
    lamprime(1)+1i*lamprime(2),vprime(1)+1i*vprime(2),para);
F3_Re = real(F3);
F3_Im = imag(F3);
F4_Re = real(F4);
F4_Im = imag(F4);

F = [F1;F2;F3_Re;F3_Im;F4_Re;F4_Im];
matlabFunction(F,'Vars',{[X;Xprime;cpseed;omega;v;lamprime;vprime];para},...
    'Optimize',false,'file','Data\F_realified.m');
DF = jacobian(F,[X;Xprime;cpseed;omega;v;lamprime;vprime]);
matlabFunction(DF,'Vars',{[X;Xprime;cpseed;omega;v;lamprime;vprime];para},...
    'Optimize',false,'file','Data\DF_realified.m');

clear F DF X Xprime rho omega v lamprime vprime para
X = sym('X',[12,1],'real');         %Symbolic variable
X0 = sym('X0',[12,1],'real');       %Taylor expansion point
DELTA = sym('DELTA',[12,1],'real'); %Arbitrary element of Delta-ball
PARA = sym('PARA',[3,1],'real');    %Symbolic parameter
DIM = 2;

F = F_realified(X,PARA);
DF = jacobian(F,X);
DF_Taylor = sym(zeros(6*DIM,6*DIM));
DF_Taylor_Rem = sym(zeros(6*DIM,6*DIM));
if ORDER == 0
    Z2_TAYLOR = subs(DF,X,X0+DELTA)-subs(DF,X,X0);
    Z2_TAYLOR_REM = sym(0);
else
    for k=1:6*DIM
        DF_Taylor(:,k) = taylor(DF(:,k),X,'order',ORDER+1,'ExpansionPoint',X0);
        DF_Taylor_Rem(:,k) = taylor(DF(:,k),X,'order',ORDER+2,'ExpansionPoint',X0) - DF_Taylor(:,k);
    end
    % Get DF(X0+DELTA)-DF(X0) and multiply through by symbolic factorial to
    % avoid problems with non-representable reciprocal floats. Get the
    % remainder and multiply by factorials too.
    Z2_TAYLOR = sym(factorial(ORDER))*(subs(DF_Taylor,X,X0+DELTA) - subs(DF,X,X0));
    Z2_TAYLOR_REM = sym(factorial(ORDER+1))*subs(DF_Taylor_Rem,X,X0+DELTA);
end
matlabFunction(Z2_TAYLOR,'Vars',{X0;DELTA;PARA},...
    'Optimize',false,'file','Data\F_Z2.m');
matlabFunction(Z2_TAYLOR_REM,'Vars',{X0;DELTA;PARA},...
    'Optimize',false,'file','Data\F_Z2_REM.m');
end