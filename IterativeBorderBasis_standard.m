function [Basis,C_0] = IterativeBorderBasis_standard(Polynomials,Variables,R,Degree,tol)

% Input: Polynomials as ususal, Variables as usual, R is represented as a
% martrix and its column vectors form the basis of R and an integer Degree

% Output: The Border Basis and also the constructed C_0 space

% this function computes the border basis iteratively

% Now we constructing C_n, starting by C_0

Number_of_Polynomials = length(Polynomials);
Monomials_of_degree_N = Construct_Monomials_standard(Degree,Variables,zeros(Variables,1)');

C_0 = [];

for i = 1:Number_of_Polynomials
    C_0 = [C_0 Coefficients(Polynomials{i},Monomials_of_degree_N)];
end


% we are now checking at which point C_n = C_{n+1}

C_n = Intersection_Spaces(Prolongation_plus(C_0,Monomials_of_degree_N,Variables),R,tol);
while Same_Spaces(C_0,C_n,tol) == 0      
     C_0 = C_n;
     C_n = Intersection_Spaces(Prolongation_plus(C_0,Monomials_of_degree_N,Variables),R,tol); %problem due to not zero elements is the intersection not always right.
end



% we are now going to contsruct the basis B we start with the 1 and extend
% it in every dicretion and check if it is needed

B = [Coefficients({{[1],[zeros(Variables,1)']}},Monomials_of_degree_N)];

B = Construct_Basis_standart(B,[Coefficients({{[1],[zeros(Variables,1)']}},Monomials_of_degree_N)],C_0,R,Variables,Monomials_of_degree_N,tol);

% What is left, is to check if B^+ is a subset of R

if rank([Prolongation_plus(B,Monomials_of_degree_N,Variables) R],tol) == rank(R,tol)
    Monomial_Basis = InverseCoefficient(B,Monomials_of_degree_N);
    Basis = [];
    for i = 1:length(Monomial_Basis)
        Basis = [Basis; Monomial_Basis{i}{1}{2}];
    end

    % avoid numerical problems

Size_C_0 = size(C_0);
for row = 1:Size_C_0(1)
    for column = 1:Size_C_0(2)
        if abs(C_0(row,column)) < tol 
            C_0(row,column) = 0;
        end
    end
end

C_0 = InverseCoefficient(C_0,Monomials_of_degree_N);

else

    % we now need to represent R in a bigger vector space, hence we need to
    % compute it anew (maybe there is a better way)

    Cardinality_of_R = size(R); % we need the column size
    Basis_Exponents = InverseCoefficient(R,Monomials_of_degree_N);
    Monomials_of_degree_larger_N = Construct_Monomials_standard(Degree + 1,Variables,zeros(Variables,1)');
    Number_of_monomials = size(Monomials_of_degree_larger_N); % number of rows is needed

    Basis_of_R = [Basis_Exponents{1}{1}{2}]; % rows represent basis vectors of R
for i = 1:Cardinality_of_R(2)
    exponent = Basis_Exponents{i}{1}{2};
    Basis_of_R = [Basis_of_R ; exponent];
    for j = 1:Variables % dividing by x_j, if possible, and adding it to R, hence making it connected to 1
        while exponent(j) - 1 >= 0
            exponent(j) = exponent(j) - 1;
            if ismember(exponent, Basis_of_R,'rows') == 0
                Basis_of_R = [Basis_of_R; exponent];
            end
        end
    end
end

    R = [];
    Cardinality_of_Basis_R = size(Basis_of_R);
    for i = 1:Cardinality_of_Basis_R(1)
    k = 1;
    while isequal(Basis_of_R(i,:),Monomials_of_degree_larger_N(k,:)) == 0
        k = k + 1;
    end
    Basis_vector = zeros(Number_of_monomials(1),1); % is already a column vector
    Basis_vector(k) = 1; % the k-th monomial in our Monomial basis Monomials_of_degree_N
    R = [R  Basis_vector];
    end

    [Basis,C_0] = IterativeBorderBasis_standard(Polynomials,Variables,Prolongation_plus(R,Monomials_of_degree_larger_N,Variables),Degree + 1,tol);

end

end