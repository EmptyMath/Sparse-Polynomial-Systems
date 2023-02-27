function Reduction_of_polynomial = Normal_Form(B,C_0,Variables,Polynomial)

% this function gets as input the monomial basis of a vector space B, as row vectors,
% with the condition 1 in B and computes for any polynomial q its normal form,
% i.e., [q] in R[x_1,...,x_n] / I. Here it is really important, that B is
% represented with a basis system!

% B is represented as row vectors

% as usual, a polynomial will be represented by
% {{[c_1],[x^[alpha_1]]},...,{[c_n],[x^[alpha_n]}}

% and as usual we also have a polynomial basis for the vector space C_0
% given in the same sense as explained above, only as a family of
% polynomials

% first, we are going to reduce every occuring monomial in q and then add
% these together

Reduction_of_monomials = {};

for i = 1:length(Polynomial)
    Reduced_monomial = Monomial_Reduction(B,C_0,Variables,Polynomial{i}{2});

    % we now multiply the polynomial by its coefficient, since we gave as
    % input above the coefficient 1

    for k = 1:length(Reduced_monomial)
        Reduced_monomial{k}{1} = Polynomial{i}{1} * Reduced_monomial{k}{1};
    end
    Reduction_of_monomials{i} = Reduced_monomial;
end

Reduction_of_polynomial = {};

% we add all reduced monomials togehter

size_of_B = size(B);

for i = 1:size_of_B(1)
    coefficient = 0;
    for k = 1:length(Reduction_of_monomials)
        l = 1;
        length_of_Reduction_monomial = length(Reduction_of_monomials{k});
        while l <= length_of_Reduction_monomial & isequal(Reduction_of_monomials{k}{l}{2},B(i,:)) == 0
            l = l +1;

        end
        if l <= length_of_Reduction_monomial
            coefficient = coefficient + Reduction_of_monomials{k}{l}{1};
        end
    end 
    Reduction_of_polynomial{i} = {[coefficient],B(i,:)};
end







end