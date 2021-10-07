function [F1,F2,F3,F4,F] = function_F(X,Xprime,cspeed,omega,v,lamprime,vprime,para)
% Implementation of the function F for the periodic traveling wave Hopf
% bifurcation. para = [rho;h;M]. The wave speed (cspeed) is the bifurcation
% parameter. Complex variable version; here, v and vprime
% are the complex unfixed coefficients of v and v', while lamprime is the
% derivative of the eigenvalue.
% para = [rho;h;M];
rho = para(1);
h = para(2);
M = para(3);
F1 = function_little_f(X,cspeed,rho);
F2 = D1f_matrix(X,cspeed,rho,h,M)*Xprime + D2f(X);
F3 = D1f_at_exp_omega_matrix(X,cspeed,rho,h,M,omega)*[1;v] - 1i*omega*[1;v];
F4 = D1f2_at_xp_and_exp_times_hv(cspeed,rho,h,M,omega,Xprime)...
    + D2D1f_at_exp_times_hv(X,cspeed,h,M,omega,v)...
    + D1f_at_theta_times_exp_omega_matrix(X,cspeed,rho,h,M,omega)*lamprime*[1;v]...
    - lamprime*[1;v] + D1f_at_exp_omega_matrix(X,cspeed,rho,h,M,omega)*[0;vprime]...
    - 1i*omega*[0;vprime];
F = [F1; F2; F3; F4];
end