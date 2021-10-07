function df = D1f2_at_exp_times_hv_matrix(c,rho,h,M,omega,xp)
df = [0; rho*xp(1)*(1-h*c/1i/omega*(1-exp(1i*omega/c))...
    -(1-h)*c/M/1i/omega*(exp(-1i*omega*M/c)-1)) ];
end