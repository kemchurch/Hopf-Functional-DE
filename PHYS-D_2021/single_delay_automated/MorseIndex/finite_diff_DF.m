function DF = finite_diff_DF(t_mesh,r,A,sigma)

h=1e-6;
N=length(t_mesh);
E=eye(N);
DF=zeros(N);

for j=1:N
    t_mesh_h = t_mesh + h*E(:,j);
    DF(:,j)=(F_mesh(t_mesh_h,r,A,sigma)-F_mesh(t_mesh,r,A,sigma))/h;
end
end