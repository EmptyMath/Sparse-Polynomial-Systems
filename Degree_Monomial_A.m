function Degree = Degree_Monomial_A(Monomial,A,Variables)

if isequal(Monomial,zeros(1,Variables)) == 1
   Degree = 0;
else
    Degree = 1;
Monomials_of_degree = Construct_Monomials_A(Degree,Variables,A);
Number_of_monomials = size(Monomials_of_degree);
while ismember(Monomial,Monomials_of_degree(1:Number_of_monomials(1),1:Variables),'rows') == 0
    Degree = Degree + 1;
    Monomials_of_degree = Construct_Monomials_A(Degree,Variables,A);
    Number_of_monomials = size(Monomials_of_degree);
end
end


end