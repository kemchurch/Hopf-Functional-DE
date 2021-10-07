function df = D1f2_at_xp_and_exp_times_hv(c,rho,h,M,omega,xp)
df = [0; rho*xp(1)*(1-h/(1i*omega)*(1-exp(1i*omega))...
    -(1-h)*(M*1i*omega)*(exp(-1i*omega*M)-1)) ];
end