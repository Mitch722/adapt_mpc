function params_set = truncate_gauss(mu, sigma, range, n)

params_set = zeros(n, length(mu));

pda11 = makedist('Normal', 'mu', mu(1), 'sigma', sigma(1));
pda22 = makedist('Normal', 'mu', mu(2), 'sigma', sigma(2));

pdb11 = makedist('Normal', 'mu', mu(3), 'sigma', sigma(3));
pdb22 = makedist('Normal', 'mu', mu(4), 'sigma', sigma(4));

for i = 1 : n
    
   params_set(i, :) = truncate
    
    
    
end