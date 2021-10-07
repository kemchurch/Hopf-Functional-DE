function [Cm,Cp,certificate] = prove_non_resonance(K1,K2,TAU,DIM,nu,N,delta1,delta2,m,meshtype,verifytype)addpath('MorseIndex')M1 = compute_M1(TAU,K1,N);M2 = compute_M2(TAU,K2,N);B = M1\eye(DIM*(N+1));M = B*M2;Eigs = eig(mid(M));if nargin==8    meshtype=0;end% C1 bound for radius 1-delta1C1m = int_compute_C1_Chebyshev(M,m,DIM,nu,1-delta1,meshtype);disp(['C1 for 1-delta1: ',num2str(C1m)]);% C1 bound for radius 1+delta2C1p = int_compute_C1_Chebyshev(M,m,DIM,nu,1+delta2,meshtype);disp(['C1 for 1+delta2: ',num2str(C1p)]);% C2 and C3 boundsif DIM==1    C2 = int_compute_C2_Chebyshev(N,B,TAU,K1,nu);    C3 = int_compute_C3_Chebyshev(M,K1,K2,TAU,nu);else    C2 = int_compute_C2_Chebyshev_v2(N,B,TAU,K1,nu);    C3 = int_compute_C3_Chebyshev_v2(M,K1,K2,TAU,nu);end% Compute C(+/-)Cm = intval(C1m)*intval(C2)*intval(C3);Cp = intval(C1p)*intval(C2)*intval(C3);if sup(Cm)<1     disp('Success in computing the index at 1-delta1.')     disp(['Product C1C2C3 = ',num2str(sup(Cm))]);     cert_m = 1;else    disp('Failure in validating the index at 1-delta1.');    disp(['Product C1C2C3 = ',num2str(sup(Cm))]);    cert_m = 0;endif sup(Cp)<1     disp('Success in computing the index at 1+delta2.')     disp(['Product C1C2C3 = ',num2str(sup(Cp))]);     cert_p = 1;else    disp('Failure in validating the index at 1+delta2.');    disp(['Product C1C2C3 = ',num2str(sup(Cp))]);    cert_p = 0;end    % Compute the indices. if strcmp(verifytype,'radpol')    % Count the number of zero eigenvalues for K2    if DIM==1        iD = K2;    else        [V,D]=eig(mid(K2));        iD = intval(zeros(DIM,1));        for k=1:DIM            [iD(k),~]=verifyeig(K2,D(k,k),V(:,k));        end    end    if sum(isnan(iD))>0        disp('Problem with the zero eigenvalues of K2.');        cert_K2_zero = 0;    else        cert_K2_zero = 1;    end    n_K2 = 0;    tol = 0;    for j=1:DIM        if ~(real(iD(j))<0 | real(iD(j))>0) & ~(imag(iD(j))<0 | imag(iD(j))<0)            n_K2 = n_K2 + 1;            tol = max(tol,sup(abs(iD(j))));        end    end    [Eigs,rad_pol] = verify_the_eigs(M,n_K2,DIM,tol+eps);    if sum(rad_pol<0)>0        cert_H = 0;        disp('Failed to verify the eigenvalues.');    else    ind_1 = count_index(Eigs,rad_pol,1-delta1);    ind_2 = count_index(Eigs,rad_pol,1+delta2);    disp(['Generalized Morse index at 1-delta1 = ',num2str(ind_1)]);    disp(['Generalized Morse index at 1+delta2 = ',num2str(ind_2)]);    if ind_1-ind_2==2        disp('Gap between radius 1-delta1 and 1+delta2 contains exactly two eigenvalues.');        cert_H = 1;    else        disp(['Error, gap between radius 1-delta1 and 1+delta2 contains ',num2str(ind_1-ind_2),' eigenvalues.']);        cert_H = 0;    end    endelseif strcmp(verifytype,'verifyeig')    cert_K2_zero = 1; %Certificate not used.    L = intval(zeros(DIM*(N+1),1));    [V,D]=eig(mid(M));    for n=1:DIM*(N+1)        [L(n),~]=verifyeig(M,D(n,n),V(:,n));    end    if sum(isnan(L)+isinf(L))>0        disp('Failed to verify the eigenvalues');        cert_H=0;    else        Ls_sup = sort(sup(abs(L)));        Ls_inf = sort(inf(abs(L)));        Lgap_sup = find(1-delta1<Ls_sup & Ls_sup<1+delta2);        Lgap_inf = find(1-delta1<Ls_inf & Ls_sup<1+delta2);        if length(Lgap_sup)==2 & length(Lgap_inf)==2            disp('Gap between radius 1-delta1 and 1+delta2 contains exactly two eigenvalues.');            cert_H=1;        else            disp('Error, gap between radius 1-delta1 and 1+delta2 contains either more than 2 eigenvalues, or some eigenvalue enclosure intersects a delta-offset circle.');            cert_H=0;        end    endendcertificate = min([cert_m,cert_p,cert_K2_zero,cert_H]);if certificate==1    disp('Success in proving non-resonance.');    % Plotting    figure    plot(real(Eigs),imag(Eigs),'*','color',[0 0 0])    hold on    color1 = [0.7 0 0];    color2 = [0.7 0.7 0];    s = (0:.002:2*pi);    circle1 = (1-delta1)*exp(1i*s); circle2 = (1+delta2)*exp(1i*s);    plot(real(circle1),imag(circle1),'color',color1,'linewidth',2)    plot(real(circle2),imag(circle2),'color',color2,'linewidth',2)    set(gca,'FontSize',20)    xlabel('$$Re(\lambda)$$','interpreter','latex','FontSize',25)    ylabel('$$Im(\lambda)$$','interpreter','latex','FontSize',25)    axis equalelse    disp('Non-resonance check failed.');end        end     