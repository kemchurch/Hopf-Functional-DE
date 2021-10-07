function C1 = int_compute_C1_Chebyshev(M,m,d,nu,r,model)

if strcmp(model,'uniform')==1
    
    theThetas = linspace(0, pi, m); 
    numTheta = length(theThetas);
    
else
     
    t_mesh = get_theta_mesh(M,r,m);
    theThetas = [0;t_mesh;pi]'; 
    numTheta = length(theThetas);
  
end
    
%%%%%%%%%%%%%%%%%%%%%%%

disp('Computing C1 ...')

Id = eye(length(M));
iM = intval(M);

C1 = 0;

for j = 1:numTheta-1
    theta = infsup(theThetas(j),theThetas(j+1));
    iA1 = inv(iM - r*exp(1i*theta)*Id);
    c1 = sup(matrix_norm(iA1,d,nu));
    C1 = max([C1,c1]);
end

C1 = max([C1,sup(1/intval(r))]);

end

