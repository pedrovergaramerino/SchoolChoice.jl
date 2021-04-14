function [test_val, pvalue, crit_val] = Test_eq_distribution(F,G,k,nperm,pctile)
% TEST_EQ_DISTRIBUTION   Use nearest-neighbor comparison to test equality of two multivariate distributions
%                        (based on W. Rolke and A. Lopez (2010) "A Test for Equality of Distributions in High Dimensions")
%                        This test is used in finding an equilibrium in Monte Carlo simulations
%
% Inputs:  
% - F: matrix of observations for first distribution
% - G: matrix of observations for second distribution
% - k: number of nearest neighbor comparisons
% - nperm: number of permutations to simulate null distribution 
% - pctile = type-I error
% 
% Output: 
% - test value = 1 if cannot reject equality of F and G at the 5% level
% - pvalue: p-value for test of equality of distributions
% - crit_val: critical value

% Sample sizes
n_F = size(F,1);
n_G = size(G,1);

% Combine multivariate samples
F_G = [F;G];

% Find k nearest neighbors for all observations in F
% N.B.: by construction, the first nearest neighbor will be the observation in F -> use k+1 neighbors to find k 
Z = knnsearch(F_G,F,'k',k+1);

% Create matrix of dummy variables: 1 if nearest neighbor is in F, 0 if in G
Z(Z>n_F)=0;
Z(Z>=1 & Z<=n_F)=1;

t_stat = sum(sum(Z(:,2:k+1)));

% Compute p-value based on binomial distribution (approximation of null distribution)
pvalue = binocdf(t_stat,n_F*k,(n_F-1)/(n_F+n_G-1),'upper');

% Find null distribution through permutations
t_stat_H0 = zeros(nperm,1);

for rr=1:nperm

    rng(rr)
    F_G_perm = F_G(randperm(size(F_G,1)),:);
    Z_perm = knnsearch(F_G_perm,F_G_perm(1:n_F,:),'k',k+1);

    Z_perm(Z_perm>n_F)=0;
    Z_perm(Z_perm>=1 & Z_perm<=n_F)=1;

    t_stat_H0(rr,1) = sum(sum(Z_perm(:,2:k+1)));

end

% Find the critical value
crit_val = prctile(t_stat_H0,pctile*100);

% return 1 if t-stat below 95th percentile, 0 otherwise
test_val = (t_stat<=crit_val);

end