function Reduced_monomial = Monomial_Reduction(B,C_0,Variables,Monomial)

% this function gets as input the monomial basis of a vector space B, as row vectors,
% with the condition 1 in B and computes for any monomial m its normal form,
% i.e., [m] in R[x_1,...,x_n] / I

% B is given as row vectors

% where C_0 is the supplementary vector space and B + C_0 = R[x] ( more or
% less), C_0 is given as {{[f_1]},...,{f_m}} and each polynomial f_i as
% {{[c_1],[x^{alpha_1}]},...,{[c_m],[x^{alpha_m}]}}

% the output is a polynomial in the form above

Index = Index_of_monomial(B,Variables,Monomial);

% if the index is 0, we are done

if Index == 0
    Reduced_monomial = {{[1],Monomial}};
    return
end

% decompose the monomial in x_i and m' such that the index of m' is stricly
% less then the index of the starting monomial

i = 0;
decompose_monomial = Monomial;

while Index_of_monomial(B,Variables,decompose_monomial) >= Index & i <= Variables
     decompose_monomial = Monomial;
     i = i + 1;
     if Monomial(i) >= 1
        decompose_monomial(i) = Monomial(i) - 1;
     end
end

reduced = Monomial_Reduction(B,C_0,Variables,decompose_monomial);

% we are now going to compute the reduction of x_i*reduced (in B^+) and then project
% it to B along C_0 

for k = 1:length(reduced)
    reduced{k}{2}(i) = reduced{k}{2}(i) + 1;
    
end

Occuring_monomials = B;
Exponents_of_C_0 = ExtractingExponents(C_0);
number_of_exponents = length(Exponents_of_C_0);

for k = 1:number_of_exponents
    if ismember(Exponents_of_C_0{k},Occuring_monomials,'rows') == 0
    Occuring_monomials = [Occuring_monomials; Exponents_of_C_0{k}];
    end
end

% constructing column basis for B and C_0

Basis_of_B = [];
cardinality_of_B = size(B);

for k = 1:cardinality_of_B(1)
    Basis_of_B = [Basis_of_B Coefficients({{[1],B(k,:)}},Occuring_monomials)];
end


Basis_of_C_0 = [];
cardinality_of_C_0 = length(C_0);

for k = 1:cardinality_of_C_0
    Basis_of_C_0 = [Basis_of_C_0 Coefficients(C_0{k},Occuring_monomials)];
end


% we are now going to project x_i*reduced onto B along C_0
x = linsolve([Basis_of_B Basis_of_C_0], Coefficients(reduced,Occuring_monomials));

projection = [];

% only the first B columns are interesting for us since its the projection 

for k = 1:cardinality_of_B(1)
    projection(k) = x(k);
end


R = InverseCoefficient(projection',B);
Reduced_monomial = R{1}; % we write R{1} because we want
% to treat the reduced monomial as a single polynomial and
% not as a family of polynomials consisting of 1 polynomial





end