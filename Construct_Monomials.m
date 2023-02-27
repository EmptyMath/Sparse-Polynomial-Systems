function Monomials_of_Degree_N = Construct_Monomials(Degree,Variables,Monomials)

% this function constructs all possible monomials of degree N in number of
% Variables and with positive or negative exponents

number_of_Monomials = size(Monomials);


% we now multiply each monomial by x_i and x_i^{-1} and check if we have
% constructed each monomial, if not we use a recursive function

for k = 1:number_of_Monomials(1)
    for i = 1:Variables
         New_Monomial_plus = Monomials(k,:);
         New_Monomial_minus = New_Monomial_plus;
         New_Monomial_plus(i) = New_Monomial_plus(i) + 1;
         if norm(New_Monomial_plus,1) <= Degree & ismember(New_Monomial_plus,Monomials,'rows') == 0
            Monomials = [Monomials ; New_Monomial_plus];
         end
         New_Monomial_minus(i) = New_Monomial_minus(i) - 1;
         if norm(New_Monomial_minus,1) <= Degree & ismember(New_Monomial_minus,Monomials,'rows') == 0
            Monomials = [Monomials ; New_Monomial_minus];
         end
    end
end

if number_of_Monomials ~= length(Monomials)
    Monomials_of_Degree_N = Construct_Monomials(Degree,Variables,Monomials);
else
    Monomials_of_Degree_N = Monomials;
end

end