function f = function_little_f(X,c,rho)
f = [X(2) ; c*X(2) - rho*X(1)*(1-X(1))];
end

