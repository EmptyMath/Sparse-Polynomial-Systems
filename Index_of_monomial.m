function Index = Index_of_monomial(Monomial_Basis,Variables,monomial)

% this function computes the index of a monomial m with respect to a
% vector space B given by a monomial basis and the condition 1 in B

% the monomial basis is given by row vectors

Index = inf;

size_1 = size(Monomial_Basis);

for n = 1:size_1(1)
    Index = min(Index,norm(( Monomial_Basis(n,:) - monomial),1));
end





end