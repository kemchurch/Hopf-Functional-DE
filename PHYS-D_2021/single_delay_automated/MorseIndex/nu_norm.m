function norm_a = nu_norm(a,nu)

N = length(a)-1;
n = (0:N)';
norm_a = sum(abs(a).*nu.^n);

end

