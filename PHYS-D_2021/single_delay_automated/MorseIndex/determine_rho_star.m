function rho_star = determine_rho_star(a,b,y1,y2)

P = [1 0 a^2*y1^2+b^2*y2^2 0 (a^2*b^2-1)*y1^2*y2^2]; 

R = roots(P);

rho_star=max(real(R));

end

