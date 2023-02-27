function Polynomial_in_Rx = Monomial_Multiplication(Polynomials,Variables)
  
% this function multiplies each polynomial by a monomial such that the new
% polynomial is contained in R[x_1,...,x_n]

% the polynomials are given in the form Polynomials {f_1,...,f_m}, and each
% polynomial is given as f_i = {{[c_1],[x^alpha_1]},...,{[c_k],[x^alpha_k]}}

Number_of_Polynomials = length(Polynomials);

for k = 1:Number_of_Polynomials
    monomial = zeros(1,Variables);
    monomials_in_polynomial_k = length(Polynomials{k});
    for l = 1:monomials_in_polynomial_k
        for i = 1:Variables
            monomial(i) = min(monomial(i),Polynomials{k}{l}{2}(i));
        end
    end
    monomial = -monomial;
    for l = 1:monomials_in_polynomial_k
        Polynomials{k}{l}{2} = Polynomials{k}{l}{2} + monomial;
    end
end

Polynomial_in_Rx = Polynomials;

end