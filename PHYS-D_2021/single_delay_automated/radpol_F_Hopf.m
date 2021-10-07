function [Y0,Z0,Z2,radpol_int,roots_double] = radpol_F_Hopf(F_HOPF,DF_HOPF,F_HOPF_Z2,F_HOPF_Z2_REM,DIM,ORDER,X,Re_FIXVAL,Im_FIXVAL,PARA,TAU,RSTAR,nrm)
Re_FIX_int = intval(Re_FIXVAL);
Im_FIX_int = intval(Im_FIXVAL);
X_int = intval(X);
PARA_int = intval(PARA);
TAU_int = intval(TAU);
RSTAR_int = intval(RSTAR);
DELTA = infsup(-1,1)*ones(6*DIM,1)*RSTAR_int;
DF = DF_HOPF(X_int,Re_FIX_int,Im_FIX_int,PARA_int,TAU_int);
A = inv(DF);
I = intval(eye(6*DIM,6*DIM));
Y0 = norm(A*F_HOPF(X_int,Re_FIX_int,Im_FIX_int,PARA_int,TAU_int),nrm);
Z0 = norm(I-A*DF,nrm);
if ORDER==-1
    Z2_main = norm(A*F_HOPF_Z2(X_int,DELTA,Re_FIX_int,Im_FIX_int,PARA_int,TAU_int),nrm);
    Z2_rem = 0;
else
    Z2_main = norm(A*F_HOPF_Z2(X_int,DELTA,Re_FIX_int,Im_FIX_int,PARA_int,TAU_int),nrm)/intval(factorial(ORDER));
    Z2_rem = norm(A*F_HOPF_Z2_REM(X_int,DELTA,Re_FIX_int,Im_FIX_int,PARA_int,TAU_int),nrm)/intval(factorial(ORDER+1));
end
Z2 = Z2_main + Z2_rem;
radpol_int = @(r)Z2*r -(1-Z0)*r + Y0;
radpol_double = [sup(Z2)-(1-sup(Z0)),sup(Y0)];
roots_double = sort(roots(radpol_double));
end