function Polynomial_Basis = InverseCoefficient(Basis,Monomials)

% this function computes a polynomial basis when a basis, where the
% basis vectors are column vectors, with coefficient vector is
% known

% the output is a family of poylnomials given by {{[f_1]},...,{[f_m]}} and
% each polynomial is again represented as a cell [f_1] =
% {{[c_1],[x^{alpha_1}]},...,{[c_m],[x^{alpha_m}]}}

size_of_Basis = size(Basis);
Cardinality_Basis = size_of_Basis(2);
Number_of_monomials = size(Monomials);
Polynomial_Basis = {};



for i = 1:Cardinality_Basis
    Number_of_monomials_in_polynomial_j = 1;
    Polynomial_vector = Basis(:,i);
    Polynomial = {};
    for k = 1:Number_of_monomials(1)
        if k <= length(Polynomial_vector) & Polynomial_vector(k) ~= 0
            Polynomial{Number_of_monomials_in_polynomial_j} = {[Polynomial_vector(k)],[Monomials(k,:)]};
            Number_of_monomials_in_polynomial_j = Number_of_monomials_in_polynomial_j + 1;
        end
    end
    Polynomial_Basis{i} = Polynomial;
    
    
end

end