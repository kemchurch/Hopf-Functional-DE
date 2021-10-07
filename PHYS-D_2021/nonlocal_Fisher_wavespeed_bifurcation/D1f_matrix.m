function dfMat = D1f_matrix(X,c,rho,h,M)
dfMat = [0,1;-rho+2*rho*X(1),c];
end