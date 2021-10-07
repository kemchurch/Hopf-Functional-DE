function sigma0 = determine_sigma0(rho_star,a,b,y1,y2)

A_star = a*y1+b*y2;
C_star = y1*y2;

sigma0 = (asin(A_star*rho_star/C_star))/2;

end

