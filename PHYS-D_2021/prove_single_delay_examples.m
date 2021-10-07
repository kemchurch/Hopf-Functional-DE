function [r_iso,Cm,Cp] = prove_single_delay_examples(NUMBER)
% This function proves the single-delay examples of Theorem [2] and Theorem
% [3]. Simply run the function with the specified input NUMBER, and the
% proof will complete. Some of them take longer than others. The output
% r_iso is the radius from the radii polynomial, and Cm and Cp are
% respectively C(1-delta1) and C(1+delta2).
BUILD=1; %always build the map F and associated data
fldr = cd;
cd('single_delay_automated'); %Navigate to the single-delay code folder
addpath('Data')
if NUMBER==1 %Lasota-Wazewska, alpha=0.3, tau=19.2089
    f = @(x,xtau,alpha,para)alpha*(-para(1)*x + exp(-xtau));
    X=[1.104542018324 0 19.208854104207  -2.703005650033   0.007171184698  -0.017941603965]';
    PARA=0.3;
    TAU=1;
    DIM=1;
    DIM_PARA=1;
    ORDER=1;
    RSTAR=1E-5;
    FIXCOMP=1;
    Re_FIXVAL=1;
    Im_FIXVAL=0;
    % Hopf isolation and transversality check
    [r_iso,cert_iso] = prove_Hopf_isolation(f,X,PARA,TAU,DIM,DIM_PARA,ORDER,RSTAR,...
        FIXCOMP,Re_FIXVAL,Im_FIXVAL,BUILD);
    if cert_iso==0
        disp('Proof failed at isolation stage.')
        cd(fldr); Cm = NaN; Cp=NaN;
        return
    end
    % Non-resonance check
    K1 = midrad(Jacobian_0(X(1),X(1),X(3),PARA),r_iso);
    K2 = midrad(Jacobian_tau(X(1),X(1),X(3),PARA),r_iso);
    delta=0.1; %+/- offset for the Morse index
    m = 20; % The mesh subdivision size of the circle in the computation of C1
    N=100; % Number of Chebyshev modes
    nu = 1.15; % The geometric decay rate for the proof
    [Cm,Cp,cert_res] = prove_non_resonance(K1,K2,DIM,TAU,nu,N,delta,delta,m,'adaptive','radpol');
    if min(cert_iso,cert_res)==1
        disp('Success, there is a Hopf bifurcation.');
    else
        disp('Hopf bifurcation proof failed.');
    end
elseif NUMBER==2 %Lasota-Wazewska, alpha=0.35, tau=37.03017
    f = @(x,xtau,alpha,para)alpha*(-para(1)*x + exp(-xtau));
    X=[1.025065556445  0  37.030171112739  -2.919994153135   0.001131897009  -0.005411625903]';
    PARA=0.35;
    TAU=1;
    DIM=1;
    DIM_PARA=1;
    ORDER=1;
    RSTAR=1E-5;
    FIXCOMP=1;
    Re_FIXVAL=1;
    Im_FIXVAL=0;
    % Hopf isolation and transversality check
    [r_iso,cert_iso] = prove_Hopf_isolation(f,X,PARA,TAU,DIM,DIM_PARA,ORDER,RSTAR,...
        FIXCOMP,Re_FIXVAL,Im_FIXVAL,BUILD);
    if cert_iso==0
        disp('Proof failed at isolation stage.')
        cd(fldr); Cm = NaN; Cp=NaN;
        return
    end
    % Non-resonance check
    K1 = midrad(Jacobian_0(X(1),X(1),X(3),PARA),r_iso);
    K2 = midrad(Jacobian_tau(X(1),X(1),X(3),PARA),r_iso);
    delta1=0.094; %+/- offset for the Morse index, lower
    delta2=0.3; %+/- offset for the Morse index, upper
    m = 25; % The mesh subdivision size of the circle in the computation of C1
    N=490; % Number of Chebyshev modes
    nu = 1.06; % The geometric decay rate for the proof
    [Cm,Cp,cert_res] = prove_non_resonance(K1,K2,DIM,TAU,nu,N,delta1,delta2,m,'uniform','radpol');
    if min(cert_iso,cert_res)==1
        disp('Success, there is a Hopf bifurcation.');
    else
        disp('Hopf bifurcation proof failed.');
    end
elseif NUMBER==3 %Coupled Lasota-Wazewska, sigma1=0.1, sigma2=0.5, tau=17, xi=0.178, unidirectional
    f=@(x,xtau,alpha,para)[-para(1)*x(1)+exp(-xtau(1)) - alpha*x(1) ; -para(2)*x(2) + exp(-xtau(2)) + alpha*x(1)];
    X = [1.14466886800860   1.08412201968491  -1.919211635684326   0.957842795546628...
         0.1780972893166  -0.154904637991545   0.868760285342032  -0.041065942913387...
         -0.047501861871352  -0.073744885404130   4.238789580398924  -6.020231755626770]';
    PARA = [0.1;0.5];
    TAU = 17;
    DIM = 2;
    DIM_PARA = 2;
    ORDER = 1;
    RSTAR=1E-5;
    FIXCOMP=1;
    Re_FIXVAL=1;
    Im_FIXVAL=0;
    % Hopf isolation and transversality check
    [r_iso,cert_iso] = prove_Hopf_isolation(f,X,PARA,TAU,DIM,DIM_PARA,ORDER,RSTAR,...
        FIXCOMP,Re_FIXVAL,Im_FIXVAL,BUILD);
    if cert_iso==0
        disp('Proof failed at isolation stage.')
        cd(fldr); Cm = NaN; Cp=NaN;
        return
    end
    % Non-resonance check
    K1 = midrad(Jacobian_0(X(1:2),X(1:2),X(5),PARA),r_iso);
    K2 = midrad(Jacobian_tau(X(1:2),X(1:2),X(5),PARA),r_iso);
    delta1=0.2; %+/- offset for the Morse index, lower
    delta2=0.2; %+/- offset for the Morse index, upper
    m = 20; % The mesh subdivision size of the circle in the computation of C1
    N=280; % Number of Chebyshev modes
    nu = 1.12; % The geometric decay rate for the proof
    [Cm,Cp,cert_res] = prove_non_resonance(K1,K2,TAU,DIM,nu,N,delta1,delta2,m,'adaptive','verifyeig');
    if min(cert_iso,cert_res)==1
        disp('Success, there is a Hopf bifurcation.');
    else
        disp('Hopf bifurcation proof failed.');
    end
elseif NUMBER==4 %Coupled Lasota-Wazewska, sigma1=0.1, sigma2=0.5, tau=17, xi=0.0677, bidirectional
    f=@(x,xtau,alpha,para)[-para(1)*x(1)+exp(-xtau(1)) - alpha*x(1) + alpha*x(2); -para(2)*x(2) + exp(-xtau(2)) + alpha*x(1) - alpha*x(2)];
    X = [1.5858236138052   0.9030535509805  -1.727579568874282   0.581589081774491...
        0.06766070440528  -0.142261377213083   0.208874995011447  -0.099748188056958...
        -0.008458308451281  -0.149652361211147   2.225857970381196  -2.762306284945626]';
    PARA = [0.1;0.5];
    TAU = 17;
    DIM = 2;
    DIM_PARA = 2;
    ORDER = 1;
    RSTAR=1E-5;
    FIXCOMP=1;
    Re_FIXVAL=1;
    Im_FIXVAL=0;
    % Hopf isolation and transversality check
    [r_iso,cert_iso] = prove_Hopf_isolation(f,X,PARA,TAU,DIM,DIM_PARA,ORDER,RSTAR,...
        FIXCOMP,Re_FIXVAL,Im_FIXVAL,BUILD);
    if cert_iso==0
        disp('Proof failed at isolation stage.')
        cd(fldr); Cm = NaN; Cp=NaN;
        return
    end
    % Non-resonance check
    K1 = midrad(Jacobian_0(X(1:2),X(1:2),X(5),PARA),r_iso);
    K2 = midrad(Jacobian_tau(X(1:2),X(1:2),X(5),PARA),r_iso);
    delta1=0.15; %+/- offset for the Morse index, lower
    delta2=0.5; %+/- offset for the Morse index, upper
    m = 20; % The mesh subdivision size of the circle in the computation of C1
    N=150; % Number of Chebyshev modes
    nu = 1.15; % The geometric decay rate for the proof
    [Cm,Cp,cert_res] = prove_non_resonance(K1,K2,TAU,DIM,nu,N,delta1,delta2,m,'adaptive','verifyeig');
    if min(cert_iso,cert_res)==1
        disp('Success, there is a Hopf bifurcation.');
    else
        disp('Hopf bifurcation proof failed.');
    end
end
cd(fldr); %Return to previous folder
end