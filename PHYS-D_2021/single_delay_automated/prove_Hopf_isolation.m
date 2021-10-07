function [r,cert] = prove_Hopf_isolation(f,X0,PARA,TAU,DIM,DIM_PARA,...
    ORDER,RSTAR,FIXCOMP,Re_FIXVAL,Im_FIXVAL,BUILD)
% [r,cert] =
% prove_Hopf_isolation(f,X0,PARA,TAU,DIM,DIM_PARA,...
%                            ORDER,RSTAR,FIXCOMP,Re_FIXVAL,Im_FIXVAL,BUILD)
% Computes the isolation radius (r) and a proof certificate (cert) for
% verification of the Hopf isolation and transversality conditions of
% Theorem 6. 
% INPUTS:
% f : anonymous function @(x,xtau,alpha,para); see documentation of
% F_Hopf_build for details.
% X0 : approximate zero of F in realified variables.
% PARA : system parameter vector; should be input as an interval enclosure
% of your actual system parameter vector. See documentation of
% F_Hopf_build for details of specification.
% TAU : system delay; should be input as an interval enclosure
% of your actual system parameters.
% DIM : positive integer, dimension of the delay DE.
% DIM_PARA : length of the PARA vector.
% ORDER : nonnegative integer, order to use for the radii polynomial.
% RSTAR : a priori maximum radius for the radii polynomial.
% FIXCOMP : component of the eigenvector that is to be fixed.
% Re_FIXVAL : real part of the fixed eigenvector component.
% Im_FIXVAL : imaginary part of the fixed eigenvector component.
% BUILD : set to 0 ONLY IF your most recent proof attempt was for the SAME
% VECTOR FIELD f, and the SAME ORDER PARAMETER was used in that proof.
% Otherwise, set to 1.
errflag = 0;
cert = 0;
if nargin==11
    BUILD = input('Build F, DF, Taylor and remainder? 1=YES, 0=NO : ');
end
if BUILD==1
    F_Hopf_build(f,DIM,DIM_PARA,FIXCOMP,ORDER);
end
[~,~,~,radpol_int,root_double] = radpol_F_Hopf(@F_Hopf,@DF_Hopf,@F_Hopf_Z2,@F_Hopf_Z2_Rem,DIM,ORDER,X0,Re_FIXVAL,Im_FIXVAL,PARA,TAU,RSTAR,inf);
if root_double<0
    disp('Proof failed! The root is nonpositive.')
    errflag = 1;
end
r = root_double+eps;
if r<intval(RSTAR)
    if radpol_int(r)<0
        disp('Radii polynomial is negative for the Hopf isolation.')
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
% Check the eigenvector is nonzero
if max(inf(abs(intval(Re_FIXVAL))),inf(abs(intval(Im_FIXVAL))))<=0
    disp('Proof failed! The eigenvector might be zero.');
    errflag = 1;
    return
end
% Check eigenvalue transversality
Re_LAM = intval(X0(4*DIM+1)) + infsup(-1,1)*r;
chk_LAM = (Re_LAM)<0 + (Re_LAM>0);
if chk_LAM==0
    disp('Proof failed! The eigenvalue crossing might be non-transversal');
    errflag = 1;
    return
end
if errflag==0
   disp('Proof of eigenvalue transversaliity successful.');
   disp(['The computed existence/uniqueness radius is: ',num2str(r),'.']);
   cert = 1;
end

end