function M2 = compute_M2(tau,K2,N)

d = length(K2);

M2 = intval(zeros(d*(N+1)));
tau = intval(tau); K2 = intval(K2);

n = intval(2:N);
i1 = intval(1); i2 = intval(2);

if d==1
    
    row_cheb = [i1 -i1/i2 -i2*(-1).^n./(n.^2-i1)];
    h_at_1 = [i1 i2*ones(1,N)];
    
    M2(1,:) = h_at_1 + (tau/i2)*K2*row_cheb;
    
    T0 = diag(intval(1:N).^(-1)*(tau/4)*K2,-1);
    T2 = -diag(intval(1:N).^(-1)*(tau/4)*K2,2);
    
    M2(2:end,:) = T0(2:end,:) + T2(1:N,1:N+1);
    
else
    
    Id = intval(eye(d));
    
    M2(1:d,1:d) = Id + (tau/2)*K2;
    M2(1:d,d+1:2*d) =  2*Id - (tau/4)*K2;
    
    for n = 2:N
        M2(1:d,n*d+1:(n+1)*d) = 2*Id - tau*(-1)^n/(n^2-1)*K2;
    end
    
    for n = 1:N-1
        M2(n*d+1:(n+1)*d,(n-1)*d+1:n*d) = (tau/(4*n))*K2;
        M2(n*d+1:(n+1)*d,(n+1)*d+1:(n+2)*d) = -(tau/(4*n))*K2;
    end
    
    M2(N*d+1:(N+1)*d,(N-1)*d+1:N*d) = (tau/(4*N))*K2;
    
end

end