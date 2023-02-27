function Exponents = ExtractingExponents(Polynomials)

% polynomials f are repesentend by a list of lists
% {{[c_alpha],[alpha]},...,{[c_beta],[beta]}} so that f{k}{1} is the
% coefficent of the k-th exponent and f{k}{2} is the k-th exponent

% Output: A list of monomials, i.e.,
% {monomial_1,...,monomial_m}, so that Exponents{k} = monomial_k
% which occurs in the polynomails given

Exponents = {};
m = length(Polynomials); % number of polynomials
k = 0;
for i = 1:m
    for l = 1:length(Polynomials{i})
        Exponents{k + l} = Polynomials{i}{l}{2};
    end
    k = k + length(Polynomials{i});
end
end