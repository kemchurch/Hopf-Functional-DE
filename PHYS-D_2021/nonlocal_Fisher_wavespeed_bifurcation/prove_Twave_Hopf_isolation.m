function r = prove_Twave_Hopf_isolation(X0,PARA,ORDER,RSTAR,BUILD)
% r = zero of the radii polynomial.
% Prove the algebraic conditions for the Hopf  bifurcation of periodic 
% traveling wave. Inputs X0 12-dimensional real, PARA 3-dimensional real,
% ORDER a natural number, RSTAR positive real.
addpath('Data')
if nargin==4
    BUILD = input('Build F, DF, Taylor and remainder? Required if this has not yet been done for specified ORDER. 1=YES, 0=NO : ');
end
DIM=2;
errflag = 0;
if BUILD==1
    disp('Building F map and related objects.');
    F_TWave_build(ORDER);
end
[~,~,~,radpol_int,root_double] = radpol_Twave(X0,PARA,RSTAR,inf,ORDER);
if root_double<0
    disp('Proof failed! The root is nonpositive.')
    errflag = 1;
end
r = root_double+eps;
if r<intval(RSTAR)
    if radpol_int(r)<0
        disp('The radii polynomial is negative at a positive r_0. Success.')
    else
        disp('Proof failed! The radii polynomial is positive.')
        errflag = 1;
        return
    end
else
    disp('Proof failed! The root is larger than RSTAR.');
    errflag = 1;
    return
end
% Check eigenvalue transversality
Re_LAM = intval(X0(4*DIM+1)) + infsup(-1,1)*r;
chk_LAM = (Re_LAM<0) + (Re_LAM>0);
if chk_LAM==0
    disp('Proof failed! The eigenvalue crossing might be non-transversal');
    errflag = 1;
    return
end
if errflag==0
   disp('Proof successful!');
   disp('The computed existence/uniqueness radius is: ');
   disp(r); 
end
end