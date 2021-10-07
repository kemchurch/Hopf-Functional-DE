function F = F_mesh(t_mesh,r,A,sigma)

m = length(t_mesh) + 2;

t_mesh = [0;t_mesh;pi];

F = zeros(m-2,1);

for j=1:m-2
    F(j) = ((h_r(r,t_mesh(j+1),sigma)+h_r(r,t_mesh(j),sigma))/2)*(t_mesh(j+1)-t_mesh(j)) - A/(m-1);
end

end

