function degree = degreePolynomial_A(Polynomial,A,Variables)

N = length(Polynomial);
degree = - inf;
for i = 1:N
    degree = max(degree,Degree_Monomial_A(Polynomial{i}{2},A,Variables));
end

end