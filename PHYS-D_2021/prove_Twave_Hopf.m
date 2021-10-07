ifunction [r_isolation,WN] = prove_Twave_Hopf(j,k,intplot)
% This function proves Wj with wave speed c_k of Theorem 4.
% If intplot = 1 (or left unspecified), the winding number plot will be
% done with intervals. Otherwise, if intplot = 0, the midpoint curve will 
% be plotted (but the proof will still be done rigorously with intervals).
% The output r_isolationis the isolation radius from
% the radii polynomial for the Hopf candidate zero
% isolation and eigenvalue transversality proof. The output WN is the
% winding number of t|->g(s(t)) for "s" the typical parameterization of the
% contour Gamma.
fldr = cd; 
cd('nonlocal_Fisher_wavespeed_bifurcation');
if j==1
    if k==1
        load('Data/symmetric_p5_2.mat');
    elseif k==2
        load('Data/symmetric_p5_1.mat');
    end
elseif j==2
    if k==1
        load('Data/symmetric_m5_2.mat');
    elseif k==2
        load('Data/symmetric_m5_1.mat');
    end
elseif j==3
    if k==1
        load('Data/leftright_047_2.mat');
    elseif k==2
        load('Data/leftright_047_1.mat');
    end
elseif j==4
    if k==1
        load('Data/GWA_8_1.mat');
    elseif k==2
        load('Data/GWA_8_2.mat');
    end
end
r_isolation = prove_Twave_Hopf_isolation(X,para,1,1E-5,1);
clear para
cd(fldr);
cd('winding_number');
ic = midrad(c,r_isolation); delta = '0.09';   epsilon = '0.05';
MESH = 3E-3;
RAY_ANGLE = 1;
if nargin<3
    intplot = 1;
end
WN = compute_winding_number(ic,h,M,rho,delta,epsilon,RAY_ANGLE,MESH,intplot);
cd(fldr);
end