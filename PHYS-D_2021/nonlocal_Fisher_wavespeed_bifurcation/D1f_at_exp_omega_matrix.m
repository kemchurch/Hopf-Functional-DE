function dfMat = D1f_at_exp_omega_matrix(X,c,rho,h,M,omega)
dfMat = [0,1 ; -rho*(1-X(1))+rho*X(1)*(-h/(1i*omega)*(1-exp(1i*omega))...
    -(1-h)/(1i*omega)/M*(exp(-1i*omega*M)-1)), c]; 
end