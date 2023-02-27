function M_lambda = Moment_Matrix_A(lambda,order,monomials,Variables,A)

% Input: a linear functional lambda, represented as a vector, with the
% columns indexed given by the monomial set monomials, represented as row
% vectors

% Output: the moment matrix of the given order corresponding to lambda

M_lambda = [];

number_of_monomials = size(monomials);

for i = 1:number_of_monomials(1)
    for j = 1:number_of_monomials(1)
        if monomials(i,Variables + 1) <= order & monomials(j,Variables + 1) <= order
        monomial = monomials(i,1:Variables) + monomials(j,1:Variables);
        
        % if Degree_Monomial_A(monomial,A,Variables) <= 2*order
            k = 1;
        while isequal(monomial,monomials(k,1:Variables)) == 0
            k = k + 1;
        end
        M_lambda(i,j) = lambda(k);
        % end
        end
    end
end




end