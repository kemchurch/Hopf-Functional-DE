function C3 = int_compute_C3_Chebyshev_v2(M,K1,K2,tau,nu)

d = length(K1);
N = length(M)/d-1;

RHO = intval(zeros(N+2,1));

tau = intval(tau); nu = intval(nu); iN = intval(N);

for n = 0:N-1
   RHO(n+1) = (1/nu^n)*tau*norm(K1*M(N*d+1:(N+1)*d,n*d+1:(n+1)*d),inf)*nu^(N+1)/(4*(iN+1)); 
end

RHO(N+1) = tau*norm(K2+K1*M(N*d+1:(N+1)*d,N*d+1:(N+1)*d),inf)*nu/(4*(iN+1)); 

RHO(N+2) = (1/nu^(iN+1))*(2+tau*norm(K2,inf)/((iN+1)^2-1)) + (1/nu)*tau*norm(K2,inf)/(4*(iN-1)) + nu*tau*norm(K2,inf)/(4*(iN+1)) ;

C3 = max(sup(RHO));

end

