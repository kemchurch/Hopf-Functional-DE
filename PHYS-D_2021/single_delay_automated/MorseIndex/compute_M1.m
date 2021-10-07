function M1 = compute_M1(tau,K1,N)

d = length(K1);

M1 = intval(zeros(d*(N+1)));
tau = intval(tau); K1 = intval(K1);

Id = intval(eye(d));

M1(1:d,1:d) = Id - (tau/2)*K1;
M1(1:d,d+1:2*d) =  (tau/4)*K1;

for n = 2:N
    M1(1:d,n*d+1:(n+1)*d) = tau*(-1)^n/(n^2-1)*K1;
end

if d==1
    
    T0 = - diag(intval(1:N).^(-1)*(tau/4)*K1,-1);
    T1 = intval(eye(N+1));
    T2 = diag(intval(1:N).^(-1)*(tau/4)*K1,2);
    M1(2:end,:) = T0(2:end,:)+T1(2:end,:)+T2(1:N,1:N+1);
    
else
    
    for n = 1:N-1
        M1(n*d+1:(n+1)*d,(n-1)*d+1:n*d) = -(tau/(4*n))*K1;
        M1(n*d+1:(n+1)*d,n*d+1:(n+1)*d) = Id;
        M1(n*d+1:(n+1)*d,(n+1)*d+1:(n+2)*d) = (tau/(4*n))*K1;
    end
 
        M1(N*d+1:(N+1)*d,(N-1)*d+1:N*d) = -(tau/(4*N))*K1;
        M1(N*d+1:(N+1)*d,N*d+1:(N+1)*d) = Id;
    
end

end