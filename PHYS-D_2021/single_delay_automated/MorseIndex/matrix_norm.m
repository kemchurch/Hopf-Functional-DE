function matrix_norm = matrix_norm(M,d,nu)

N = length(M)/d-1;

if d == 1
    tM = M;
    nu_power = nu.^abs(0:N);
    matrix_norm = max(sum(abs(tM).*repmat(nu_power,N+1,1)')./nu_power);
else
    
    % We define tLambda as an (N+1)x(N+1) matrix containting the ||*||_d norm of each block
    tM = intval(zeros(N+1));
    
    for m = 0:N
        for n = 0:N
            tM(m+1,n+1) = norm(M(m*d+1:(m+1)*d,n*d+1:(n+1)*d),inf);
        end
    end
    
    nu_power = nu.^abs(0:N);
    matrix_norm = max(sum(abs(tM).*repmat(nu_power,N+1,1)')./nu_power);
    
end

end