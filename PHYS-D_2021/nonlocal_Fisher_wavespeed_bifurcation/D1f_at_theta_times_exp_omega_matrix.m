function dfMat = D1f_at_theta_times_exp_omega_matrix(X,c,rho,h,M,omega)
bigterm = h*(1+exp(1i*omega)*(1i*omega-1))/omega^2 + ...
    (1-h)*exp(-1i*omega*M)*((1-exp(1i*omega*M))+1i*omega*M)/M/omega^2;
dfMat = [0,1 ; -rho*(1-X(1))+rho*X(1)*bigterm, c];
end