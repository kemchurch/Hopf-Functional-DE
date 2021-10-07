function []=ZeroPlots_cspeed(rho,h,M,c_range,omega_range)
% para = [rho;h;M];
f = @(c,omega)det(D1f_at_exp_omega_matrix([1;0],c,rho,h,M,omega) - 1i*omega*eye(2,2));
fimplicit(@(rho,omega)real(f(rho,omega)),[c_range(1),c_range(2),omega_range(1),omega_range(2)],'r');
hold on
fimplicit(@(rho,omega)imag(f(rho,omega)),[c_range(1),c_range(2),omega_range(1),omega_range(2)],'b');
end