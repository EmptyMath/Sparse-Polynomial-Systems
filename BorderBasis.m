function [Basis,C_0] = BorderBasis(Polynomials,variables,tol)

% polynomials f are repesentend by a list of lists
% {{[c_alpha],[alpha]},...,{[c_beta],[beta]}} so that f{k}{1} is the
% coefficent of the k-th exponent and f{k}{2} is the k-th exponent

% make sure that every monomial is contained in R[x_1,...,x_n]

% variables is the number of variables used

% extracting all exponents given occurring in the polynomials in question
Exponents = ExtractingExponents(Polynomials);


N = maxdegree(Exponents); %maybe not important
Number_of_monomials = size(Construct_Monomials(N + 1,variables,zeros(variables,1)')); % number of monomials of degree at most N

% constructing the vector space R mentioned in Algorithm 2,
% we will construct a basis which is connected to 1, by considering 
% each exponent contained in Exponents and dividing it by x_j if possible

M = length(Exponents);
Basis_of_R = [Exponents{1}]; % rows represent basis vectors of R
for i = 1:M
    exponent = Exponents{i};
    if ismember(exponent, Basis_of_R,'rows') == 0
        Basis_of_R = [Basis_of_R; exponent];
    end
    for j = 1:variables % dividing by x_j, if possible, and adding it to R, hence making it connected to 1
        while exponent(j) - 1 >= 0
            exponent(j) = exponent(j) - 1;
            if ismember(exponent, Basis_of_R,'rows') == 0
                Basis_of_R = [Basis_of_R; exponent];
            end
        end
    end
end


% Constructing R and C_0 as vector spaces, therefore we map every exponent
% to a specfic vector 

Number_of_basis_elements = length(Basis_of_R);
R = [];
Monomials_of_degree_N = Construct_Monomials(N + 1,variables,zeros(variables,1)'); % we have to consider N+1 since we will also multiply some polynomials by x_i and x_i^{-1}

% Firstly R, now every column vector will be a basis vector:
% we need to take into account that also negative exponents can occur!!!


for i = 1:Number_of_basis_elements
    k = 1;
    while isequal(Basis_of_R(i,:),Monomials_of_degree_N(k,:)) == 0
        k = k + 1;
    end
    Basis_vector = zeros(Number_of_monomials(1),1); % is already a column vector
    Basis_vector(k) = 1; % the k-th monomial in our Monomial basis Monomials_of_degree_N
    R = [R  Basis_vector];
end


[Basis,C_0] = IterativeBorderBasis(Polynomials,variables,R,N+1,tol);

end