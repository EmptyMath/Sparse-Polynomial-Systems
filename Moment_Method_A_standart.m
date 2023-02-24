function Real_Roots = Moment_Method_A_standard(Polynomials,A,Variables,tol)

% this function computes the real roots of a given polynomial system, where
% the polynomials are given in the form Polynomials {f_1,...,f_m}, and each
% polynomial is given as f_i = {{[c_1],[x^alpha_1]},...,{[c_k],[x^alpha_k]}}
% Variables is an integer

% This algorithm only works for A subset of N^n

% tol discribes a tolerance for the rank computation

% step (2) + (3)

Exponents = ExtractingExponents(Polynomials);
D = - inf;
K = length(Polynomials);
for k = 1:K
    D = max(D,degreePolynomial_A(Polynomials{k},A,Variables));
end

t = D

T = 0;
Monomials = Construct_Monomials_A(t,Variables,A);

while T == 0
    t
    Monomials = Construct_Monomials_A(t,Variables,A);
    [At,c,b,K] = SDP_Input_A(Polynomials,Monomials,A,Variables,t);
    

    pars.eps = 1e-4;
    pars.beta = 0.1;
    pars.theta = 0.1;
    [x,~,info] = sedumi(At,b,c,K,pars)
    info
    
    Number_of_Monomials = size(Monomials)
    size_c = size(c)
    %mat(x(Number_of_Monomials(1) + 1:size_c(1)))
    %rank(mat(x(Number_of_Monomials(1) + 1:size_c(1))),tol)
    [Condition,order] = Rank_Condition_A(D,floor(t/2),x,Monomials,tol,Variables,A)
    if Condition == 0 
    t = t + 1;
    else
        % step (4)

        order
        M_order = Moment_Matrix_A(x,order,Monomials,Variables,A)
        rank(M_order,tol)
        size(M_order)
        

Size_M = size(M_order);
rank(M_order,tol)
if rank(M_order,tol) == Size_M(1)
    Basis_J = zeros(Size_M(1),1);
else
    Basis_J = null(M_order,tol);
end

% avoid nuermical problems

Size_J = size(Basis_J);
for column = 1:Size_J(2)
    for row = 1:Size_J(1)
        if abs(Basis_J(row,column)) < tol
            Basis_J(row,column) = 0;
        end
    end
end


% we compute a Basis of the quotient space corresponding to J and also the
% multiplication operators

% First, we need to represent the Basis of J as polynomials and not as
% vectors

Number_of_Monomials_now = size(Monomials);
Monomials_Now = Monomials(1:Number_of_Monomials_now(1),1:Variables);
Monomials = Monomials_Now;

Polynomial_Basis = InverseCoefficient(Basis_J,Monomials);

% we now make use of the classical Border Basis, hence we need to construct
% new polynomials




% we have to make sure that every polynomial in the polynomial basis is
% contained in R[x_1,...,x_n]. Therefore we need to multiply each
% polynomial by some monomial.


Polynomial_Basis = Monomial_Multiplication(Polynomial_Basis,Variables);

% We are now computing a border basis in order to compute the
% multiplication maps and apply the eigenvalue method

[Border_Basis,C_0] = BorderBasis_standart(Polynomial_Basis,Variables,tol);
Size_of_Border_Basis = size(Border_Basis)
Size_Border = size(Border_Basis);
if Size_Border(1) == rank(M_order,tol) & info.numerr < 2 & rank(At) == length(c) % avoid numerical problems
   T = 1;
else
    t = t + 1;
end
    end
end

% avoid nuermical problems

Size_C_0 = size(C_0);
for row = 1:Size_C_0(1)
    for column = 1:Size_C_0(2)
        if abs(C_0(row,column)) < tol 
            C_0(row,column) = 0;
        end
    end
end

C_0_Polynomials = InverseCoefficient(C_0,Monomials);

Multiplication_Maps = Multiplication_Operators(Border_Basis,C_0_Polynomials,Variables)

% Eigenvalue method
Multiplication_Random = zeros(Size_of_Border_Basis(1),Size_of_Border_Basis(1));

for i = 1:Variables
    Random_Coefficient = rand;
    Multiplication_Random = Multiplication_Random + Random_Coefficient * Multiplication_Maps(1:Size_of_Border_Basis(1),(i-1) * Size_of_Border_Basis(1) + 1:i * Size_of_Border_Basis(1))
end

Real_Roots = Eigenvalue_Method(Border_Basis,Multiplication_Random,Multiplication_Maps,Variables)
end