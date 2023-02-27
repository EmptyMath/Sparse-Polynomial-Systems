function B_plus = Prolongation_plus(B, Monomials,Variables)

% constructes the prolongation of a basis B, where the basis vectors are
% column vectors

Cardinality_Basis = size(B);

% we use InverseCoefficient in order to multiply these polynomials by x_i
% and then will translate them back to vectors

Monomial_Basis_B = InverseCoefficient(B,Monomials);


for i = 1:Variables
    for k = 1:Cardinality_Basis(2)
        New_polynomial_vector = Monomial_Basis_B{k};
        Number_of_Monomials = length(New_polynomial_vector);
        
        % multiply the polynomial by x_i
        New_Monomial_vector_plus = New_polynomial_vector;
       for l = 1:Number_of_Monomials
           New_Monomial_vector_plus{l}{2}(i) = New_Monomial_vector_plus{l}{2}(i) + 1;
       end
       
       Coefficient_plus = Coefficients(New_Monomial_vector_plus,Monomials);
       

        if ismember(Coefficient_plus',B','rows') == 0
            B = [B Coefficient_plus];
        end

        
    end 

end

% we now get a generator system of B^+ (but is does not need to be a basis)
B_plus = B;


end
