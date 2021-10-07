function r=enc_eig(A,V,lambda)

n = size(V,1);
iV = intval(V);

[m,IV] = max(abs(V));

O = ones(n,1);
df = intval(A)-intval(lambda)*speye(n);

F = df*iV;

df(:,IV) = -iV;

R = inv(mid(df));
Y = sup(abs(R*F));
Z0 = sup(abs(eye(n)-R*df))*O;

O(IV) = 0;
Z1 = 2*abs(R)*O;

% coefficients of radii polynomials [r^2,r,const]
rd = [Z1,Z0-ones(n,1),Y];

%INT=evaluate_intersection_int(rd); % evaluate where all the polynomials are negative
r = evaluate_intersection_int(rd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
 function neg_inter=evaluate_intersection_int(rd)
            
            neg_inter=NaN;
            Rm=(-rd(:,2)-sqrt(rd(:,2).^2-4*rd(:,1).*rd(:,3)))./(2*rd(:,1));
            Rp=-Rm-rd(:,2)./rd(:,1);
            ISCOMP=find(isreal(Rp)+isreal(Rm)<2, 1);
            if isempty(ISCOMP)==0%isempty(rd(ISCOMP,1)>1)==0
                disp('error: a polynomial is positive everywhere')
                r = -inf;
                return
            end
            
            intersection=[max(Rm),min(Rp)];
            if intersection(1)<intersection(2)
                neg_inter=intersection(1);
            end
                                                       
 end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   