function [Y0,Z0,Z2,radpol_int,roots_double] = radpol_Twave(X,PARA,RSTAR,nrm,ORDER)
if nargin==2
    RSTAR=1E-4;
    nrm=inf;
    ORDER=2;
end
addpath('Data')
DIM = 2;
X_int=intval(X);
PARA_int = intval(PARA);
RSTAR_int = intval(RSTAR);
DELTA = infsup(-1,1)*ones(6*DIM,1)*RSTAR_int;
DF = DF_realified(X_int,PARA_int);
A = inv(DF);
I = intval(eye(6*DIM,6*DIM));
Y0 = norm(A*F_realified(X_int,PARA_int),nrm);
Z0 = norm(I-A*DF,nrm);
Z2_main = norm(A*F_Z2(X_int,DELTA,PARA_int),nrm)/intval(factorial(ORDER));
Z2_rem = norm(A*F_Z2_REM(X_int,DELTA,PARA_int),nrm)/intval(factorial(ORDER+1));
Z2 = Z2_main + Z2_rem;
radpol_int = @(r)(Z2+Z0-1)*r + Y0;
radpol_double = [-1+sup(Z2)+sup(Z0),sup(Y0)];
roots_double = sort(roots(radpol_double));
end