function C2 = int_compute_C2_Chebyshev_v2(N,iB,tau,K1,nu)

d = length(K1);

v1 = intval(zeros(N+1,1)); v2 = intval(zeros(N+1,1));

for m = 0:N
    v1(m+1) = norm(tau*(-1^(N+1))/((N+1)^2-1)*iB(m*d+1:(m+1)*d,1:d)*K1+tau/(4*N)*iB(m*d+1:(m+1)*d,N*d+1:(N+1)*d)*K1,inf);
    v2(m+1) = norm(iB(m*d+1:(m+1)*d,1:d)*K1,inf);
end

nu = intval(nu); tau = intval(tau); K1 = intval(K1); iN = intval(N);

rho_N = tau*norm(K1,inf)*nu/(4*(iN+1)) ;
rho_Nplus1 = (1/nu^(iN+1))*nu_norm(v1,nu) + tau*norm(K1,inf)*nu/(4*(iN+2)) ;
rho_infty = (1/nu^(iN+2))*tau/((iN+2)^2-1)*nu_norm(v2,nu) + (1/nu)*tau*norm(K1,inf)/(4*(iN+1)) + nu*tau*norm(K1,inf)/(4*(iN+3)) ;

rho = intval(max(sup([rho_N rho_Nplus1 rho_infty])));

C2 = sup(max([sup(matrix_norm(iB,d,nu)),1])/(1-rho));

end

