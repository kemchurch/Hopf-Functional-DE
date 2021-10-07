function C2 = int_compute_C2_Chebyshev(N,iB,tau,K1,nu)

nu = intval(nu); tau = intval(tau); K1 = intval(K1); iN = intval(N);
 
rho_N = tau*abs(K1)*nu/(4*(iN+1)); 

v = (tau*K1*(-1)^(iN+1)/((iN+1)^2-1))*iB(:,1) + tau*K1/(4*iN)*iB(:,N+1);

rho_Nplus1 = (1/nu^(iN+1))*nu_norm(v,nu) + tau*abs(K1)*nu/(4*(iN+2)); 

rho_infty = (1/nu^(iN+2))*tau*abs(K1)/((iN+2)^2-1)*nu_norm(iB(:,1),nu) + (1/nu)*tau*abs(K1)/(4*(iN+1)) + nu*tau*abs(K1)/(4*(iN+3)) ;

rho = intval(max(sup([rho_N rho_Nplus1 rho_infty])));

C2 = sup(max([sup(matrix_norm(iB,1,nu)),1])/(1-rho));

end

