function ind = get_mesh(C1_tot,numTheta)

V = abs(C1_tot(2:end)-C1_tot(1:end-1)); % Variation of C1
N = numTheta;
I = sum(V); % Integral of V over [0,pi]
figure, plot((1:N)*pi/N,V,'linewidth',3), axis tight

m = 170; % Desired size of the adaptive mesh in theta

ind(1) = 1;

for j = 2:m+100
    k=0;
    while (ind(j-1)+k)< N && (sum(V(ind(j-1):ind(j-1)+k)) < I/m)
        k=k+1;
        if (ind(j-1)+k) == N
            ind(j) = numTheta;
            return
        end
    end
    ind(j) = ind(j-1) + k;
end

end

