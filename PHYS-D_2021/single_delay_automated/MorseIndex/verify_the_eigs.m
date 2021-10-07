function [Eigs,rad_pol]=verify_the_eigs(iA,m,d,tolerance)

%%% This code takes an interval matrix (possibly of radius 0), computes its
%%% eigendecomposition using the "eig" function in MATLAB, and then
%%% validates the existence of each eigenpair (v,lambda) of the nonlinear
%%% problem iA*v-lambda*v=0. It does so using the radii polynomials. If
%%% successful, the code returns the list of eigenvalues Eigs together
%%% with the radius enclosing them. Hence, the interval
%%% [Eigs(k)-rad_pol,Eigs(k)+rad_pol] constains a unique eigenvalue

A = mid(iA);

[V,D] = eig(A);
Eigs = diag(D); 

index = find(abs(Eigs)>=tolerance);

Eigs = Eigs(index);
V = V(:,index);
n = length(index);

N = length(A)/d - 1;

rad_pol = zeros(1,n);

if n + m ~= d*(N+1)
    disp('problem with the zeros eigenvalues of K2')
    return
end

for j = 1:n
    rad_pol(j)=enc_eig(iA,V(:,j),Eigs(j));
end

end