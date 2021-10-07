clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           CHOOSE A MODEL          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% model = 1 : Mackey-Glass          %%%
%%% model = 2 : Cubic Ikeda-Matsumoto %%%
%%% model = 3 : Delayed van der Pol   %%%
%%% model = 4 : Predator-Prey model   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = 4;

%%%%%%%%%%%%%%%%%%%%
%%% Mackey-Glass %%%
%%%%%%%%%%%%%%%%%%%%

if model == 1
    
    tau = 2; gamma = 1; beta = 2; rho = 10; K1 = -gamma;
    %c = 0; nu = 1.3; N = 32; d = 1; r = 1; m = 10;
    %c = 0; nu = 1.2; N = 70; d = 1; r = .6; m = 10;
    %c = 0; nu = 1.15; N = 200; d = 1; r = .29; m = 40;
    %c = 0; nu = 1.05; N = 600; d = 1; r = .2; m = 60;
    %c = 1; nu = 1.1; N = 310; d = 1; r = 1; m = 100; % can't find an adaptive mesh
    c = 1; nu = 1.2; N = 130; d = 1; r = .85; m = 30;
    %c = 1; nu = 1.1; N = 250; d = 1; r = .46; m = 30;
    %c = 1; nu = 1.05; N = 500; d = 1; r = .341; m = 90;
    K2 = beta*(1+(1-rho)*c^rho)/((1+c^rho)^2); n = 0; % No zero eigenvalues

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cubic Ikeda-Matsumoto %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if model == 2
    
    tau = 1.59; K1 = 0; n = 0; % No zero eigenvalues in K2
    %c = 0; nu = 1.6; N = 5; K2 = 1-3*c^2; d = 1; r = 1; m = 10;
    %c = 0; nu = 1.2; N = 50; K2 = 1-3*c^2; d = 1; r = .25; m = 20;
    %c = 1; nu = 1.9; N = 6; K2 = 1-3*c^2; d = 1; r = 1; m = 20;
    c = 1; nu = 1.25; N = 72; K2 = 1-3*c^2; d = 1; r = .31; m = 25;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Delayed van der Pol %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if model == 3
    
    % c = (0,0) is the only fixed point
    tau = 2; kappa = -1; epsilon = 0.15;
    K1 = [[0 1];[kappa epsilon]];
    K2 = [[0 0];[-1 0]];
    N = 14; nu = 1.5; c = [0;0]; d = 2; n = N; r = 1; m = 25; % N repeated zero eigenvalues
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Predator-prey model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if model == 4
    
    %(y1,y2) is the unique positive equilibrium
    r1 = 2; r2 = 1; a = 1; b = 1/2; 
    y1 = (r2+b*r1)/(a*b+1); y2 = (r1-a*r2)/(a*b+1);
    rho_star = determine_rho_star(a,b,y1,y2);
    sigma0 = determine_sigma0(rho_star,a,b,y1,y2);
    tau = sigma0/rho_star-.1;
    K1 = [[-a*tau*y1 0];[0 -b*tau*y2]];
    K2 = [[0 -tau*y1];[tau*y2 0]];
    N = 25; nu = 1.25; c = [0;0]; d = 2; r = 1; n = 0; m = 20; % 0 repeated zero eigenvalues
    % can't find an adaptive mesh
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% c = 0; tau = 1; N = 5; nu = 2.4; r = 1;
% d = 12; ind = 2:d-1; m = 15; %works for d = 3,6,12
% %d = 24; ind = 2:d-1; m = 25; 
% K1 = diag(log([3 2 1./ind])); K2 = 1e-2*rand(d); 

M1 = compute_M1(tau,K1,N);
M2 = compute_M2(tau,K2,N);

B = M1\intval(eye(d*(N+1))); % Computation of the exact inverse of M1

M = B*M2; % M^N = (M_1^N)^(-1)*M_2^N

%%%%%%%%%%%%%
%%%   C1  %%%
%%%%%%%%%%%%%

C1 = int_compute_C1_Chebyshev(M,m,d,nu,c,r,model);
disp(['C1 = ',num2str(C1)])

%%%%%%%%%%%%%
%%%   C2  %%%
%%%%%%%%%%%%%

d = length(K1);

if d==1
    C2 = int_compute_C2_Chebyshev(N,B,tau,K1,nu);
else
    C2 = int_compute_C2_Chebyshev_v2(N,B,tau,K1,nu);
end

disp(['C2 = ',num2str(C2)])

%%%%%%%%%%%%%
%%%   C3  %%%
%%%%%%%%%%%%%

if d==1
    C3 = int_compute_C3_Chebyshev(M,K1,K2,tau,nu);    
else
    C3 = int_compute_C3_Chebyshev_v2(M,K1,K2,tau,nu);
end

disp(['C3 = ',num2str(C3)])

%%%%%%%%%%%%
%%%   C  %%%
%%%%%%%%%%%%

C = intval(C1)*intval(C2)*intval(C3);

disp(['C = ',num2str(sup(C))])

if sup(C)<1
    disp('Success in computing the index')
else
    disp('Failure in computing the index')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Spectrum of M  %%%
%%%%%%%%%%%%%%%%%%%%%%%%

disp('Rigorously computing the spectrum of M ...')
[Eigs,rad_pol] = verify_the_eigs(M,n,d);
index = count_index(Eigs,rad_pol,r);
disp(['Generalized Morse index = ',num2str(index)])

figure
plot(real(Eigs),imag(Eigs),'*','color',[0 0 0])
hold on
color = [0.7 0 0];
s = (0:.002:2*pi); circle = r*exp(1i*s); plot(real(circle),imag(circle),'color',color,'linewidth',2)

set(gca,'FontSize',20)

xlabel('$$Re(\lambda)$$','interpreter','latex','FontSize',25)
ylabel('$$Im(\lambda)$$','interpreter','latex','FontSize',25)

axis equal