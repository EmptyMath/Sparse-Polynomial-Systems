function Monomials_of_Degree_N = Construct_Monomials_A(Degree,Variables,A)

% this function constructs all possible monomials of degree N by the set
% A, given by row vectors where the last entry is the degree of the given
% monomial,i.e., [Monomial Deg_A(Monomial)]

number_of_Monomials = size(A);
Monomials_of_Degree_N = [];
for k = 1:number_of_Monomials(1)
    if isequal(A(k,:),zeros(1,Variables)) == 1
    Monomials_of_Degree_N = [Monomials_of_Degree_N ; A(k,:) 0];
    else
    Monomials_of_Degree_N = [Monomials_of_Degree_N ; A(k,:) 1];
    end

end


% we now add all elements in A as often as the degree 

if Degree >= 2
    for d = 2:Degree
        Monomials_of_Degree_N = Adding_A(Monomials_of_Degree_N,A,d,Variables);
    
    end
end



end