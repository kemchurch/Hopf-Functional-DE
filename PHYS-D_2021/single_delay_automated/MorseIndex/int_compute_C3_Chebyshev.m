function C3 = int_compute_C3_Chebyshev(M,K1,K2,tau,nu)

N = length(M)-1;

RHO = intval(zeros(N+2,1));

tau = intval(tau); nu = intval(nu); iN = intval(N);

for n = 0:N-1
   RHO(n+1) = (1/nu^n)*tau*abs(K1*M(N+1,n+1))*nu^(N+1)/(4*(iN+1)); 
end

RHO(N+1) = tau*abs(K2+K1*M(N+1,N+1))*nu/(4*(iN+1)); 

RHO(N+2) = (1/nu^(iN+1))*(2+tau*abs(K2)/((iN+1)^2-1)) + (1/nu)*tau*abs(K2)/(4*(iN-1)) + nu*tau*abs(K2)/(4*(iN+1)) ;

C3 = max(sup(RHO));

end

