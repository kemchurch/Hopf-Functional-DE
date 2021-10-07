function [t_mesh] = get_theta_mesh(M,r,m)

sigma = eig(mid(M)); sigma = sigma(1:min(20,length(M)));

%%% Fine computation of the integral of h_r over [0,pi]
N = 1000; 
t_mesh_M = linspace(0, pi, N);
vec = zeros(1,N-1);

for k = 1:N-1
    vec(k) = h_r(r,t_mesh_M(k),sigma);
end

Int = pi/(N-1)*sum(vec);

%%% Computation of the mesh of size m on which each integrals is equal

tol = 1e-13; it = 0;

t_mesh = linspace(0, pi, m);
t_mesh = t_mesh(2:m-1)';

F = F_mesh(t_mesh,r,Int,sigma);

%display(['||F|| = ', num2str(norm(F,inf))])

while norm(F)>tol && it<20

    DF = finite_diff_DF(t_mesh,r,Int,sigma);
    t_mesh = t_mesh - DF\F_mesh(t_mesh,r,Int,sigma);
    F=F_mesh(t_mesh,r,Int,sigma);
    it = it + 1;
    %display(['||F|| = ', num2str(norm(F,inf))])
    
end

if norm(F)>tol || min(t_mesh(2:end)-t_mesh(1:end-1))<0
    disp('Failure with the adaptive mesh computation')
else
%    disp('Success with the adaptive mesh computation')
end

end

