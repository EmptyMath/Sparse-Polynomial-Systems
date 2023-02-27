function Degrees = Degree_of_Polynomials_A(Polynomials,Variables,A)

N = length(Polynomials);

Degrees = [];
for i = 1:N
    Degrees = [Degrees degreePolynomial_A(Polynomials{i},A,Variables)];
end

end