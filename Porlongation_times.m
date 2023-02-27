function C_n_times = Porlongation_times(Basis_of_C_0, Monomials,Variables)

% this function constructes the porlongation of a vector space by the
% coordinates of x_i and x_i^{-1}

Cardinality_Basis = length(Basis_of_C_0);
Polynomial_Basis_C_0 = InverseCoefficient(Basis_of_C_0,Monomials);
for i = 1:Variables
    for k = 1:Cardinality_Basis
        New_polynomial_vector = Polynomial_Basis_C_0{k};
        Monomials_in_polynomial = length(New_polynomial_vector);
        
        % multiply the polynomial by x_i
        New_polynomial_vector_plus = New_polynomial_vector;

       for l = 1:Monomials_in_polynomial
           New_polynomial_vector_plus{l}{2}(i) = New_polynomial_vector_plus{1}{2}(i) + 1;
       end
       
       Coefficient_plus = Coefficients(New_polynomial_vector_plus,Monomials);

        if ismember(Coefficient_plus',Basis_of_C_0','rows') == 0
            Basis_of_C_0 = [Basis_of_C_0 Coefficient_plus];
        end

         % multiply the polynomial by x_i^{-1}
        New_polynomial_vector_minus = New_polynomial_vector;

       for l = 1:Monomials_in_polynomial
           New_polynomial_vector_minus{l}{2}(i) = New_polynomial_vector_minus{1}{2}(i) - 1;
       end
       
       Coefficient_minus = Coefficients(New_polynomial_vector_minus,Monomials);

        if ismember(Coefficient_minus',Basis_of_C_0','rows') == 0
            Basis_of_C_0 = [Basis_of_C_0 Coefficient_minus];
        end
    end 

end

C_n_times = Basis_of_C_0;


end