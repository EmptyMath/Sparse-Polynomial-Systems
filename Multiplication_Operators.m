function Matrices = Multiplication_Operators(Basis,C_0,Variables);

% this function computes the multiplication operators m_i for i =
% 1,...,Variables in the form of a single matrix so M = [m_1,...,m_n]

% C_0 is represented as a family of polynomials and the Basis as row
% vectors

Matrices = [];
Size_of_Basis = size(Basis);



for i = 1:Variables
    
    polynomial = zeros(1,Variables);
    polynomial(i) = 1;

    for k = 1:Size_of_Basis(1)
       Matrices = [Matrices Coefficients(Normal_Form(Basis,C_0,Variables,{{[1],Basis(k,:) + polynomial}}),Basis)];
end

end