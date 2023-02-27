function Coefficient_Vector = Coefficients(Polynomial,Monomials)

% this function computes the coefficient vector of a polynomial f given in
% the form of {{[c_alpha],[alpha]},...,{[c_beta],[beta]}} so that f{k}{1} is the
% coefficent of the k-th exponent and f{k}{2} is the k-th exponent with
% regards to given Monomials (check that each Monomial in f occurs in the
% set Monomials !!!)

% the monomials are given as row vectors 


% checking if each Monomial in f occurs also in the set Monomials

Number_of_monomials = length(Polynomial);
for i = 1:Number_of_monomials
    if ismember(Polynomial{i}{2},Monomials,'rows') == 0
        Coefficient_Vector = [];
    end
end

cardinality = size(Monomials);
Coefficient_Vector = zeros(cardinality(1),1);

% we are now looking which monomials occur in the basis at which index

for i = 1:Number_of_monomials
    for k = 1:cardinality(1)
        if isequal(Polynomial{i}{2},Monomials(k,:)) == 1
            Coefficient_Vector(k) = Coefficient_Vector(k) + Polynomial{i}{1}; % or just + 1instead of + coefficient ? does not change the vector space
        end
    end
end


end