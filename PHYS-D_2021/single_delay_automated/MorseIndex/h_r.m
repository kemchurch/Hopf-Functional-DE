function h = h_r(r,theta,sigma)

h = 1/min(abs(sigma-r*exp(1i*theta)));

end

