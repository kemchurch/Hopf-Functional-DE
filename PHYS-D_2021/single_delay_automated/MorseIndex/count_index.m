function index= count_index(Eigs,rad_pol,r)

iEigs = intval(Eigs);

N = inf(abs(iEigs) + intval(rad_pol'));

index = length(N(N>r));

end

