
close all
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Mackey-Glass %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 50;
tau = 2; gamma = 1; beta = 2; rho = 10; K1 = -gamma; M1 = compute_M1(tau,K1,N);

%%%%%%%%%%%%%
%%% c = 0 %%%
%%%%%%%%%%%%%

c = 0 ; K2 = beta*(1+(1-rho)*c^rho)/((1+c^rho)^2); M2 = compute_M2(tau,K2,N); M = mid(M1\M2);
[~,D] = eig(M); Eigs = diag(D); 

figure
plot(real(Eigs),imag(Eigs),'*','color',[0 0 0])
hold on

set(gca,'FontSize',20)

xlabel('$$Re(\lambda)$$','interpreter','latex','FontSize',25)
ylabel('$$Im(\lambda)$$','interpreter','latex','FontSize',25)

axis equal
box off

r = 1; color = [0.7 0 0]; s = (0:.002:2*pi); circle = r*exp(1i*s); plot(real(circle),imag(circle),'color',color,'linewidth',2)
r = .6; color = [0 0.7 0]; s = (0:.002:2*pi); circle = r*exp(1i*s); plot(real(circle),imag(circle),'color',color,'linewidth',2)
r = .29; color = [0 0 0.7]; s = (0:.002:2*pi); circle = r*exp(1i*s); plot(real(circle),imag(circle),'color',color,'linewidth',2)
r = .2; color = [0.8500, 0.3250, 0.0980]; s = (0:.002:2*pi); circle = r*exp(1i*s); plot(real(circle),imag(circle),'color',color,'linewidth',2)

%%%%%%%%%%%%%
%%% c = 1 %%%
%%%%%%%%%%%%%

c = 1 ; K2 = beta*(1+(1-rho)*c^rho)/((1+c^rho)^2); M2 = compute_M2(tau,K2,N); M = mid(M1\M2);
[V,D] = eig(M); Eigs = diag(D); 

figure
plot(real(Eigs),imag(Eigs),'*','color',[0 0 0])
hold on

set(gca,'FontSize',20)

%title('Mackey-Glass spectrum of $$DF(1)$$','interpreter','latex','FontSize',25)
xlabel('$$Re(\lambda)$$','interpreter','latex','FontSize',25)
ylabel('$$Im(\lambda)$$','interpreter','latex','FontSize',25)

axis equal
box off

r = 1; color = [0.7 0 0]; s = (0:.002:2*pi); circle = r*exp(1i*s); plot(real(circle),imag(circle),'color',color,'linewidth',2)
r = .85; color = [0 0.7 0]; s = (0:.002:2*pi); circle = r*exp(1i*s); plot(real(circle),imag(circle),'color',color,'linewidth',2)
r = .46; color = [0 0 0.7]; s = (0:.002:2*pi); circle = r*exp(1i*s); plot(real(circle),imag(circle),'color',color,'linewidth',2)
r = .341; color = [0.8500, 0.3250, 0.0980]; s = (0:.002:2*pi); circle = r*exp(1i*s); plot(real(circle),imag(circle),'color',color,'linewidth',2)
